classdef QuadprogSolver < IQpSolver
%QUADPROGSOLVER QP Solver based on MATLAB quadprog.

properties(Access = private)
    options;    % Option for the solver
end

methods
    function setup(obj)
    % Setup the solver
        obj.options = optimoptions('quadprog',...
            'display','off', ...
            'MaxIterations',1000, ...           % Default: 200
            'OptimalityTolerance', 1e-10, ...   % Default: 1e-8
            'StepTolerance', 1e-16, ...         % Default: 1e-12
            'ConstraintTolerance', 1e-8 ...     % Default: 1e-8
        );
    end
    
    function [x,cost,exitflag] = solve(obj,H,f,A,b,lb,ub)
    % Solve a QP problem defined as follow
    % min 1/2*x'*H*x + f'*x   subject to:  A*x <= b  and lb <= x <= ub
        [x,cost,exitflag] = quadprog(...
            H,f,A,b,[],[],lb,ub,[],obj.options...
        );
    end
    
    function reset(~)
    % Reset the solver
        
    end
    
    function release(~)
    % Release the solver
        
    end
end

methods(Static)
    function errMessage = getWarningMessage(exitflag)
    % Return error message corresponding to the exitflag. If there is no
    % error, return an empty message.
        errMessage = '';
        if exitflag < 1
            switch exitflag
                case 0
                    errMessage = [errMessage, ...
                        'Maximum number of iterations exceeded.']; 
                case -2
                    errMessage = [errMessage, ...
                        'Problem is infeasible.'];
                case -3
                    errMessage = [errMessage, ...
                        'Problem is unbounded.'];
                case -6
                    errMessage = [errMessage, ...
                        'Nonconvex problem detected.'];
                case -8
                    errMessage = [errMessage, ...
                        'Unable to compute a step direction.'];
                otherwise
                    errMessage = [errMessage, 'The exitflag is ', ...
                        num2str(exitflag)];
            end
        end
    end
end

end

