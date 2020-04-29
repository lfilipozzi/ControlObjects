classdef IQpSolver < handle
% IQPSOLVER Interface for QP solvers.

methods(Abstract)
    setup(obj)
    % Setup the solver
    
    [x,cost,exitflag] = solve(obj,H,f,A,b,lb,ub)
    % Solve a QP problem defined as follow
    % min 1/2*x'*H*x + f'*x   subject to:  A*x <= b and lb <= x <= ub
    
    reset(obj)
    % Reset the solver
    
    release(obj)
    % Release the solver
    
end

methods(Static)
    function errMessage = getWarningMessage(~)
    % Return error message corresponding to the exitflag. If there is no
    % error, return an empty message.
        errMessage = '';
    end
end

end

