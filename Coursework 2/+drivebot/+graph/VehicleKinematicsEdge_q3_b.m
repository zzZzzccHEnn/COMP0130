classdef VehicleKinematicsEdge_q3_b < g2o.core.BaseBinaryEdge
    properties(Access = protected)
        dT;
        Odom
    end

    methods(Access = public)
        function this = VehicleKinematicsEdge_q3_b(dT)
            assert(dT >= 0);
            this = this@g2o.core.BaseBinaryEdge(3);            
            this.dT = dT;
        end

        function initialize(this)
            priorX = this.edgeVertices{1}.x;

            c = cos(priorX(3));
            s = sin(priorX(3));
            
            M = [c -s 0;
                s c 0;
                0 0 1];

            dx = this.edgeVertices{2}.x - this.edgeVertices{1}.x;
            dx(3) = g2o.stuff.normalize_theta(dx(3));
            this.Odom = inv(M) * (dx/this.dT);

            % Compute the posterior assming no noise
%             this.edgeVertices{2}.x = this.edgeVertices{1}.x + M * this.Odom;

            % Wrap the heading to -pi to pi
%             this.edgeVertices{2}.x(3) = g2o.stuff.normalize_theta(this.edgeVertices{2}.x(3));

        end

        function computeError(this)
    
            % Q1b:
            % Complete implementation
%             warning('vehiclekinematicsedge:computeerror:unimplemented', ...
%                     'Implement the rest of this method for Q1b.');

%             ################################

            priorX = this.edgeVertices{1}.x;

            c = cos(priorX(3));
            s = sin(priorX(3));
            
            M = [c s 0;
                -s c 0;
                0 0 1];            

            dx = this.edgeVertices{2}.x - priorX;
            dx(3) = g2o.stuff.normalize_theta(dx(3));
%             this.errorZ = M * dx / this.dT - this.z;
            this.errorZ = M * dx / this.dT - this.Odom;
            this.errorZ(3) = g2o.stuff.normalize_theta(this.errorZ(3));

%             ################################

        end
        
        % Compute the Jacobians
        function linearizeOplus(this)

            % Q1b:
            % Complete implementation
%             warning('vehiclekinematicsedge:linearizeoplus:unimplemented', ...
%                 'Implement the rest of this method for Q1b.');

%             ################################
            priorX = this.edgeVertices{1}.x;

            c = cos(priorX(3));
            s = sin(priorX(3));
            
            M = [c s 0;
                -s c 0;
                0 0 1];

            dx = this.edgeVertices{2}.x - priorX;
            this.J{2} = M/this.dT;
            this.J{1}(1, 1) = - c / this.dT;
            this.J{1}(1, 2) = - s / this.dT;
            this.J{1}(1, 3) = -dx(1) * s / this.dT + dx(2) * c / this.dT;
            this.J{1}(2, 1) = s / this.dT;
            this.J{1}(2, 2) = - c / this.dT;
            this.J{1}(2, 3) = -dx(1) * c / this.dT - dx(2) * s / this.dT;
            this.J{1}(3, 3) = -1 / this.dT;

%             ################################

        end
        


    end
end