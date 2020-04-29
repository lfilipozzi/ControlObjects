classdef QpOasesSolver < IQpSolver
%QPOASESSOLVER qpOASES solver.

properties(Access = private)
    options;    % Option for the solver
end

methods
    function setup(obj)
    % Setup the solver
        obj.options = qpOASES_options('reliable');
    end
    
    function [x,cost,exitflag] = solve(obj,H,f,A,b,lb,ub)
    % Solve a QP problem defined as follow
    % min 1/2*x'*H*x + f'*x   subject to:  A*x <= b  and lb <= x <= ub
        [x,cost,exitflag] = qpOASES(H,f,A,lb,ub,[],b,obj.options);
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
        if exitflag ~= 0
            switch exitflag
                case 1
                    errMessage = [errMessage, ...
                        'Maximum number of iterations exceeded.']; 
                case -1
                    errMessage = [errMessage, ...
                        'Internal error.'];
                case -2
                    errMessage = [errMessage, ...
                        'Problem is unfeasible.'];
                case -3
                    errMessage = [errMessage, ...
                        'Problem is unbounded.'];
                otherwise
                    errMessage = [errMessage, 'The exitflag is ', ...
                        num2str(exitflag)];
            end
        end
    end
end

end

