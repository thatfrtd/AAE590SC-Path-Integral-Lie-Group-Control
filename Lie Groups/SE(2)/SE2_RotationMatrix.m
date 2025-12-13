classdef SE2_RotationMatrix < Group
    %SE2_ROTATIONMATRIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name = "SE2 Rotation Matrix"
        size = 9
        dim = 3
        element
        identity = eye(3)
    end
    
    methods
        function obj = SE2_RotationMatrix(R, r)
            arguments
                R = eye(2)
                r = zeros([2, 1])
            end
            %SE2_ROTATIONMATRIX Construct an instance of this class
            %   Detailed explanation goes here
            obj.element = [R, r; zeros([1, 2]), 1];
        end

        function val = R(obj)
            val = obj.element(1:2, 1:2);
        end
        function val = r(obj)
            val = obj.element(1:2, 3);
        end
        
        function composition = compose(X, Y)
            composition = SE2_RotationMatrix(X.R * Y.R, X.r + X.R * Y.r);
        end
        function inverse = inv(X)
            inverse = X;
            inverse.element = [X.R', -X.R' * X.r;
                               0, 0, 1];
        end

        function val = constraint(X)
            val = X.R' * X.R;
        end
        function val = act(X, p)
            val = X.r + X.R * p;
        end

        function tau = vee(tau_hat)
            theta = tau_hat(2, 1);
            rho = tau_hat(1:2, 3);
            tau = [rho; theta];
        end
        function tau_hat = hat(G, tau)
            rho = tau(3);
            theta = tau(1:2);
            theta_cross = skew2(theta);
            tau_hat = G;
            tau_hat.element = [theta_cross, rho; 
                               0, 0, 0];
        end

        function X = Exp(G, tau)
            theta = tau(1);
            rho = tau(2:3);
            R = make_R2(theta);
            r = G.V(theta) * rho;
            X = SE2_RotationMatrix(R, r);
        end
        function tau = Log(X)
            R = X.R;
            r = X.r;
            theta = atan2(R(2, 1), R(1, 1));
            rho = X.V_inv(theta) * r;
            tau = [rho; theta];
        end

        function obj = cayley(obj, X)

        end
        function obj = inv_cayley(obj, X)

        end

        function Adjoint = Ad(X)
            % Manifold adjoint
            Adjoint = [X.R, -skew2(1) * X.r;
                       0, 0, 1];
        end
        
        function val = V(obj, theta)
            % SO(2) left jacobian
            if theta ~= 0
                val = sin(theta) / theta * eye(2) + (1 - cos(theta)) / theta * skew2(1);
            else
                val = eye(2); % limit of sin(x) / x as x -> 0 is 1
            end
        end
        function val = V_inv(obj, theta)
            % SO(2) inverse left jacobian
            val = sin(theta) / theta * eye(2) + (1 - cos(theta)) / theta * skew2(1);
        end
    end

    methods (Static)
        function adjoint = ad(G, twist)
            % Tangent space adjoint
            adjoint = [zeros([1, 3]); -skew2(1) * twist(2:3), skew2(twist(1))]; % Verify??
        end
    end
end

