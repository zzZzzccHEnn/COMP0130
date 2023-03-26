% This class implements an event-based estimation system using g2o and
% the barebones for building up a minimal, ideal SLAM system. The system is
% event-based and responds to a sequence of events which are time stamped
% and served in order. To implement your SLAM system, you will need to
% implement various methods which mostly involve working with the graph.
% These methods are initially stubbed out and will generate exceptions if
% you try to call them.

classdef DriveBotSLAMSystem < minislam.slam.SLAMSystem
    
    properties(Access = public, Constant)
        % Platform state dimension
        NP = 3;
        
        % Landmark dimension
        NL = 2;
        
        % Initial cache size; might help a bit with performance
        INITIAL_CACHE_SIZE = 10000;
    end
    
    properties(Access = protected)
        
        % The most recently created vehicle vertex.
        currentVehicleVertex;
        
        % The set of all vertices associated with the vehicle state over
        % time.
        vehicleVertices;
        vehicleVertexId;
        
        % The set of all prediction edges. These are removed from the graph
        % afterwards if we don't use prediction
        processModelEdges;
        numProcessModelEdges;
        
        % The landmark vertices. Confusingly enough, "Map" here refers to
        % the data structure which is used to store the landmarks. (It
        % allows random access of landmarkID to landmark object.)
        landmarkIDStateVectorMap;
        
        % How often we recommend running the optimization
        recommendOptimizationPeriod;
        
        % Flag to show if we should prune the edges. This is needed for
        % question Q3a
        removePredictionEdgesFromGraph;
        keepFirstPredictionEdge;

        % Flag to show if we should prune the edges
        graphPruning;
        
    end
    
    methods(Access = public)
        
        % Create the localization system and start it up.
        function this = DriveBotSLAMSystem(configuration)
            
            % Call the base class constructor
            this = this@minislam.slam.SLAMSystem(configuration);
            
            % Preallocate for convenience
            this.vehicleVertices = cell(1, this.INITIAL_CACHE_SIZE);
            
            % No vehicle vertices initally set
            this.vehicleVertexId = 0;
            
            % The set of prediction edges, initially empty
            this.processModelEdges = cell(1, this.INITIAL_CACHE_SIZE);
            this.numProcessModelEdges = 0;
            
            % Allocate the landmark map
            this.landmarkIDStateVectorMap = containers.Map('KeyType', 'int64', 'ValueType', 'any');
            
            % By default, run very infrequently
            this.recommendOptimizationPeriod = inf;
            
            this.removePredictionEdgesFromGraph = false;
            this.keepFirstPredictionEdge = false;

            this.graphPruning = false;
        end
        
        % Destroy the graph when we destroy the SLAM system.
        % Without this, MATLAB will crash whenever this object is destroyed.

        function delete(this)
            vertices = this.graph.vertices();

            for v = 1 : length(vertices)
                this.graph.removeVertex(vertices{v});
            end
        end
        
        % Recommend if an optimization is a good idea. Based on an event,
        % some activities (e.g., such as loop closing) can have a very big
        % impact on the estimates. The logic we have here just recommends
        % an optimization if a fixed number of steps have been completed.
        
        function recommendation = recommendOptimization(this)
            
            % This is how to do it after every 100 steps
            recommendation = rem(this.stepNumber, ...
                this.recommendOptimizationPeriod) == 0;
        end
        
        % Set the value of how often recommend optimization should return
        % true
        function setRecommendOptimizationPeriod(this, newRecommendOptimizationPeriod)
            this.recommendOptimizationPeriod = newRecommendOptimizationPeriod;
        end
        
        % Return the current mean and covariance estimate of the robot.
        % This is only valid after optimization has been called.
        function [x, P] = platformEstimate(this)
            [xS, PS] = this.graph.computeMarginals(this.currentVehicleVertex);
            x=full(xS);
            P=full(PS);
        end
        
        % Returns the entire history of the platform estimates. Suppose
        % there are n vehicle vertices. T is a 1 by N dimensional vector of
        % timesteps. X is a 3 by N dimensional vector of vehicle state (x,
        % y, theta). P is a 3 by N dimensional vector where the nth column
        % are the diagonals from the covariance matrix.
        function [T, X, P] = platformEstimateHistory(this)
            
            % Extract the graph
            [xS, PS] = this.graph.computeMarginals();
            
            % Create the output array
            X = zeros(this.NP, this.vehicleVertexId);
            P = zeros(this.NP, this.vehicleVertexId);
            T = zeros(1, this.vehicleVertexId);
            
            % Copy the outputs over
            for v = 1 : this.vehicleVertexId
                idx = this.vehicleVertices{v}.hessianIndex();
                
                T(v) = this.vehicleVertices{v}.time();
                
                % Copy the estimate into the array. If the vertices is
                % fixed (conditioned), its estimate is okay. The covariance
                % is not explicitly defined, but has a value of zero.
                % Therefore we fill this manually.
                if (isempty(idx) == true)
                    X(:, v) = this.vehicleVertices{v}.estimate();
                    P(:, v) = zeros(3, 1);
                else
                    X(:, v) = full(xS(idx));
                    P(:, v) = full(diag(PS(idx, idx)));
                end
            end
        end
        
        % Return the means and covariances of the landmark estimates. These
        % are only valid after optimization has been called.
        function [x, P, landmarkIds] = landmarkEstimates(this)
            
            landmarkVertices = values(this.landmarkIDStateVectorMap);
            
            numberOfLandmarks = length(landmarkVertices);
            
            landmarkIds = NaN(1, numberOfLandmarks);
            x = NaN(this.NL, numberOfLandmarks);
            P = NaN(this.NL, this.NL, numberOfLandmarks);
            
            [xS, PS] = this.graph.computeMarginals();
            
            for l = 1 : numberOfLandmarks
                landmarkIds(l) = landmarkVertices{l}.landmarkId();
                idx = landmarkVertices{l}.hessianIndex();
                x(:, l) = full(xS(idx));
                if (isempty(idx == true))
                    P(:, :, l) = zeros(3, 3);
                else
                    P(:, :, l) = full(PS(idx, idx));
                end
            end
        end
        
        % We overload the optimize method so that you can add additional
        % logic here
        function chi2 = optimize(this, maximumNumberOfOptimizationSteps)
            
            % Remove the prediction edges if requested.
            if (this.removePredictionEdgesFromGraph == true)
                this.deleteVehiclePredictionEdges();
            end
            
            % Now call the actual optimizer. Let it handle the default if
            % no steps are specified.
            if (nargin > 1)
                chi2 = optimize@minislam.slam.SLAMSystem(this, ...
                    maximumNumberOfOptimizationSteps);
            else
                chi2 = optimize@minislam.slam.SLAMSystem(this);
            end
        end
        
        function setRemovePredictionEdges(this, removeEdges, keepFirst)
            this.removePredictionEdgesFromGraph = removeEdges;
            this.keepFirstPredictionEdge = keepFirst;
            
        end

        function setQ3_b(this, pruning)
            this.graphPruning = pruning;
        end
    end
    
    % These are the methods you will need to overload
    methods(Access = protected)
        
        % Handle the initial condition
        
        function handleInitialConditionEvent(this, event)
            
            % Create the first vertex, set its estimate to the initial
            % value and add it to the graph.
            this.currentVehicleVertex = drivebot.graph.VehicleStateVertex(this.currentTime);
            this.currentVehicleVertex.setEstimate(event.data);
            this.graph.addVertex(this.currentVehicleVertex);
            
            % Set the book keeping for this initial vertex.
            this.vehicleVertexId = 1;
            this.vehicleVertices{this.vehicleVertexId} = this.currentVehicleVertex;
            
            % If the covariance is 0, the vertex is known perfectly and so
            % we set it as fixed. If the covariance is non-zero, add a
            % unary initial prior condition edge instead. This adds a soft
            % constraint on where the state can be.
            if (det(event.covariance) < 1e-6)
                this.currentVehicleVertex.setFixed(true);
            else
                initialPriorEdge = drivebot.graph.InitialPriorEdge();
                initialPriorEdge.setMeasurement(event.data);
                initialPriorEdge.setInformation(inv(event.covariance));
                initialPriorEdge.setVertex(this.currentVehicleVertex);
                this.graph.addEdge(initialPriorEdge);
            end
        end
        
        function handleNoPrediction(~)
            % Nothing to do
        end
        
        function handleHeartbeatEvent(this, ~)
            % Nothing to do
        end
        
        function handlePredictToTime(this, time, dT)

            % Create the next vehicle vertex and add it to the graph
            
            this.currentVehicleVertex = drivebot.graph.VehicleStateVertex(time);
            
            % Q1b:
            % Implement prediction code here
%             warning('drivebotslam:handlepredicttotime:unimplemented', ...
%                 'Implement the rest of this method for Q1b.');

%             #############################
%             set the initial estimation to the vertex and add to the graph
          
            processEdge = drivebot.graph.VehicleKinematicsEdge(dT);
            processEdge.setVertex(1,this.vehicleVertices{this.vehicleVertexId});
            processEdge.setVertex(2,this.currentVehicleVertex);
            processEdge.setMeasurement(this.u);
            processEdge.setInformation(inv(this.uCov));
            processEdge.initialize();
            this.graph.addEdge(processEdge);

            this.currentVehicleVertex.setEstimate(this.currentVehicleVertex.estimate());
            this.graph.addVertex(this.currentVehicleVertex);

            this.numProcessModelEdges = this.numProcessModelEdges + 1;
            this.processModelEdges{this.numProcessModelEdges} = processEdge;

%             #############################

            
            % Bump the indices
            this.vehicleVertexId = this.vehicleVertexId + 1;
            this.vehicleVertices{this.vehicleVertexId} = this.currentVehicleVertex;
        end
        
        function handleGPSObservationEvent(this, event)

            % Q1c:
            % Create a GPS measurement edge and add it to the graph
%             warning('drivebotslam:handlegpsobservationevent:unimplemented', ...
%                 'Implement the rest of this method for Q1c.');
            GPSmeasurementEdge = drivebot.graph.GPSMeasurementEdge(this.configuration.gpsPositionOffset);
            GPSmeasurementEdge.setVertex(1,this.currentVehicleVertex);
            GPSmeasurementEdge.setMeasurement(event.data);
            GPSmeasurementEdge.setInformation(inv(event.covariance));
            this.graph.addEdge(GPSmeasurementEdge);


        end
        
        function handleCompassObservationEvent(this, event)
            
            % Create a compass measurement edge and add it to the graph
            compassMeasurementEdge = drivebot.graph.CompassMeasurementEdge(this.configuration.compassAngularOffset);
            compassMeasurementEdge.setVertex(1, this.currentVehicleVertex);
            compassMeasurementEdge.setMeasurement(event.data);
            compassMeasurementEdge.setInformation(inv(event.covariance));
            this.graph.addEdge(compassMeasurementEdge);
        end
        
        function handleLandmarkObservationEvent(this, event)
            
            % Iterate over all the landmark measurements
            for l = 1 : length(event.landmarkIds)
                
                % Get the landmark vertex associated with this measurement.
                % If necessary, a new landmark vertex is created and added
                % to the graph.
                [landmarkVertex, newVertexCreated] = this.createOrGetLandmark(event.landmarkIds(l));
                z = event.data(:, l);

                % Q2b:
                % Complete the implementation
%                 warning('drivebotslamsystem:handlelandmarkobservationevent:unimplemented', ...
%                     'Implement the rest of this method for Q2b.');
%             #############################
                landmarkRangeBearingEdge = drivebot.graph.LandmarkRangeBearingEdge();
                landmarkRangeBearingEdge.setVertex(1,this.currentVehicleVertex);
                landmarkRangeBearingEdge.setVertex(2,landmarkVertex);
                landmarkRangeBearingEdge.setMeasurement(z);
                landmarkRangeBearingEdge.setInformation(inv(event.covariance));

                if (newVertexCreated == true)
                    landmarkRangeBearingEdge.initialize();
                end
%             #############################
                
                this.graph.addEdge(landmarkRangeBearingEdge);
            end
        end
        
        function deleteVehiclePredictionEdges(this)

            % Q3a:            
%             warning('drivebotslam:deletevehiclepredictionedges:unimplemented', ...
%                 'Implement the rest of this method for Q3a.');
%             #############################

            if (this.graphPruning == false)
                    %             get all the edges from the graph
                allEdges = this.graph.edges();
%               initialise the number of deleted edges and the first
%               prediction edge
                deleted = 0;
                numbfirstPrediction = -1;

%               keep the first prediction edge or not and then initialise the
%               flag
                if (this.keepFirstPredictionEdge == true)
                    firstPrediction = false;
                else
                    firstPrediction = true;
                end

%               loop all the edges in the graph
                for e=1:length(allEdges)
%                   if keep the first prediction edge, store the index of the
%                   first prediction edge
                    if (class(allEdges{e}) == "drivebot.graph.VehicleKinematicsEdge")
                        if firstPrediction == false
                            firstPrediction = true;
                            numbfirstPrediction = e;
                        end
                    
%                       remove the edges
                        if (e ~= numbfirstPrediction)
                            this.graph.removeEdge(allEdges{e});
                            deleted = deleted + 1;
                        end
                    end
                end
                disp("####################################");
                disp(deleted);
            else
                Pruning(this);
            end


%             #############################

        end

        function Pruning(this)
            % Q3_b
            fprintf('the number of edges before pruning = %.2f. \n', length(this.graph.edges()));

            allVertices = this.vehicleVertices;
            % initialise the threshold 
            count = 0;
            threshold = 0.02;

            % Loop the vehivle verteces
            for v = 1:length(allVertices)-1
                if (isempty(allVertices{v})==false && isempty(allVertices{v+1})==false && isempty(allVertices{v+2})==false)
                    vertex1 = allVertices{v}.estimate();
                    vertex2 = allVertices{v+1}.estimate();
                    vertex3 = allVertices{v+2}.estimate();
                    % Calculate the angular velocity of the heading and
                    % compare with the threshold
                    % "0.1" is the dT
                    deltaHeading_1 = (vertex2(3)-vertex1(3))/0.1;
                    deltaHeading_2 = (vertex3(3)-vertex2(3))/0.1;
                    if (abs(deltaHeading_1) < threshold && abs(deltaHeading_2) < threshold)
                        count = count + 1;
                        thisVertex = allVertices{v+1};
                        edgesWithThisVertex = thisVertex.edges();
                        % Remove all the landmark edges related to this
                        % vertex
                        for e = 1:length(edgesWithThisVertex)
                            if (class(edgesWithThisVertex{e})=="drivebot.graph.LandmarkRangeBearingEdge")
                                this.graph.removeEdge(edgesWithThisVertex{e});
                            end
                        end
                    end
                end
            end


            % Retrieve the vertices without connecting to any landmark
            isThereLandmark = false;
            verticesWithoutLandmark = [];
            vehivleVertices = [];
            allVertices = this.graph.vertices();
            for v = 2:length(allVertices)
                if (class(allVertices{v}) == "drivebot.graph.VehicleStateVertex")
                    vehivleVertices = [vehivleVertices v];
                    edgesWithThisVertex = allVertices{v}.edges();
                    for e = 1:length(edgesWithThisVertex)
                        if (class(edgesWithThisVertex{e}) == "drivebot.graph.LandmarkRangeBearingEdge")
                            isThereLandmark = true;
                            continue
                        end
                    end
                end
                if (isThereLandmark == false)
                    verticesWithoutLandmark = [verticesWithoutLandmark v];
                end
                isThereLandmark = false;
            end

            diffVertices = setdiff(vehivleVertices, verticesWithoutLandmark);

            % delete the edges of the vehicle vertices without connecting
            % any landmark
            for v = 3:length(vehivleVertices)
                edgesWithThisVertex = allVertices{vehivleVertices(v)}.edges();
                for e = 1:length(edgesWithThisVertex)
                    if (class(edgesWithThisVertex{e}) == "drivebot.graph.VehicleKinematicsEdge")
                        this.graph.removeEdge(edgesWithThisVertex{e});
                    end
                end
            end

            diffVertices = [2 diffVertices];

            % create the edge to connect the vehicle vertices with
            % connection to the landmarks
            for i = 1:length(diffVertices)-1
                Vertex_1 = allVertices{diffVertices(i)};
                Vertex_2 = allVertices{diffVertices(i+1)};
                processEdge = drivebot.graph.VehicleKinematicsEdge_q3_b(0.1*(diffVertices(i+1)-diffVertices(i)));
                processEdge.setVertex(1,Vertex_1);
                processEdge.setVertex(2,Vertex_2);
                processEdge.setMeasurement(this.u);
                processEdge.setInformation(inv(this.uCov));
                processEdge.initialize();
                this.graph.addEdge(processEdge);
            end

%             disp(diffVertices);
           
            fprintf('the number of edges after pruning = %.2f. \n', length(this.graph.edges()));
            fprintf('the number of redundant vertecs = %.2f. \n', count);

        end
        
        
        % This method returns a landmark associated with landmarkId. If a
        % landmark exists already, it is returned. If it does not exist, a
        % vertex is created and is added to the graph.
        function [landmarkVertex, newVertexCreated] = createOrGetLandmark(this, landmarkId)

            
            % If the landmark exists already, return it
            if (isKey(this.landmarkIDStateVectorMap, landmarkId) == true)
                landmarkVertex = this.landmarkIDStateVectorMap(landmarkId);
                newVertexCreated = false;
                return
            end
            
            fprintf('Creating landmark %d\n', landmarkId);
            
            % Create the new landmark add it to the graph
            landmarkVertex = drivebot.graph.LandmarkStateVertex(landmarkId);
            this.landmarkIDStateVectorMap(landmarkId) = landmarkVertex;
            
            this.graph.addVertex(landmarkVertex);
            
            newVertexCreated = true;
        end
        
        function storeStepResults(this)
            % Nothing
        end
        
    end
end
