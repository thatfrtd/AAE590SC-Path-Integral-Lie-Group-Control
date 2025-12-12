classdef SE3_RotationMatrix < Group
    %SE3_ROTATIONMATRIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name = "SE3 Rotation Matrix"
        size = 16
        dim = 6
        element
        identity = eye(4)
    end
    
    methods
        function obj = SE3_RotationMatrix(R, r)
            arguments
                R = eye(3)
                r = zeros([3, 1])
            end
            %SE3_ROTATIONMATRIX Construct an instance of this class
            %   Detailed explanation goes here
            obj.element = [R, r; zeros([1, 3]), 1];
        end

        function val = R(obj)
            val = obj.element(1:3, 1:3);
        end
        function val = r(obj)
            val = obj.element(1:3, 4);
        end
        
        function composition = compose(X, Y)
            composition = SE3_RotationMatrix(X.R * Y.R, ...
                                             X.r + X.R * Y.r);
        end
        function inverse = inv(X)
            inverse = X;
            inverse.element = [X.R', -X.R' * X.r;
                               0, 0, 0, 1];
        end

        function val = constraint(X)
            val = X.R' * X.R;
        end
        function val = act(X, p)
            val = X.r + X.R * p;
        end

        function tau = vee(tau_hat)
            theta = [tau_hat(3, 2); tau_hat(1, 3); tau_hat(2, 1)];
            rho = tau_hat(1:3, 4);
            tau = [rho; theta];
        end
        function tau_hat = hat(G, tau)
            rho = tau(1:3);
            theta = tau(4:6);
            theta_cross = skew(theta);
            tau_hat = G;
            tau_hat.element = [theta_cross, rho; 
                               0, 0, 0, 0];
        end

        function X = Exp(G, tau)
            theta = tau(1:3);
            rho = tau(4:6);
            theta_mag = norm(theta);
            theta_dir = theta / theta_mag;
            if theta_mag ~= 0
                R = eye(3) + sin(theta_mag) * skew(theta_dir) + (1 - cos(theta_mag)) * skew(theta_dir) ^ 2;
            else
                R = eye(3);
            end
            X = SE3_RotationMatrix(R, ...
                                   G.V(theta_mag, theta_dir) * rho);
        end
        function tau = Log(X)
            theta_mag = acos((trace(X.R) - 1) / 2);
            Rdiff = X.R - X.R';
            Rdiff_vee = [Rdiff(3, 2); Rdiff(1, 3); Rdiff(2, 1)];
            theta = theta_mag * Rdiff_vee / (2 * sin(theta_mag));
            rho = X.V_inv(theta_mag, theta) * X.r;
            tau = [rho; theta];
        end

        function obj = cayley(obj, X)

        end
        function obj = inv_cayley(obj, X)

        end

        function Adjoint = Ad(obj, X)
            Adjoint = [X.R, -skew(X.r) * X.R;
                       0, 0, 0, X.R];
        end

        function val = V(obj, theta_mag, theta_dir)
            % SO(3) left jacobian
            if theta_mag ~= 0
                val = eye(3) + (1 - cos(theta_mag)) / theta_mag * skew(theta_dir) + (theta_mag - sin(theta_mag)) / theta_mag * skew(theta_dir) ^ 2;
            else
                val = eye(3); % limit of sin(x) / x as x -> 0 is 1
            end
        end
        function val = V_inv(obj, theta_mag, theta_dir)
            % SO(3) inverse left jacobian
            val = eye(3) - 1/2 * theta_mag * skew(theta_dir) + (1 - theta_mag * (1 + cos(theta_mag)) / (2 * sin(theta_mag))) * skew(theta_dir) ^ 2;
        end
    end

    methods (Static)
        function adjoint = ad(G, twist)
            % Adjoint map
            w = twist(1:3);
            v = twist(4:6);

            adjoint = [skew(w), zeros(3,3); skew(v), skew(w)];
        end
    end
end

