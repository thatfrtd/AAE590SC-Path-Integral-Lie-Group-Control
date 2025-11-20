classdef R2 < Group
    %R2 Summary of this class goes here
    %   Detailed explanation goes here
    %   Scalar last quaternions

      properties
        Name = "R2"
        size = 2
        element
        dim = 2
        identity = [0; 0]
    end
    
    methods
        function obj = R2(v)
            arguments
                v = [0; 0]
            end
            %R2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.element = v;
        end

        function val = v(obj)
            val = obj.element;
        end
        
        function composition = compose(X, Y)
            composition = R2(X.element + Y.element);
        end
        function inverse = inv(X)
            inverse = R2(-X.element);
        end

        function val = constraint(X)
            val = X.compose(X.inv);
        end
        function val = act(X, p)
            val = X.element + p;
        end

        function tau = vee(G, tau_hat)
            tau = tau_hat;
        end
        function tau_hat = hat(G, tau)
            tau_hat = tau;
        end

        function X = Exp(G, tau)
            X = R2(tau);
        end
        function tau = Log(X)
            tau = X.element;
        end

        function val = Ad(X)
            val = eye(X.dim);
        end

        function val = ad(X, twist)
            val = eye(X.dim);
        end

        function X = cayley(G, tau)
            
        end
        function tau = inv_cayley(X)
            
        end
    end
end

