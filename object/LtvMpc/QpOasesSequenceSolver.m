classdef QpOasesSequenceSolver < IQpSolver
%QPOASESSEQUENCESOLVER qpOASES_sequence solver.

properties(Access = private)
    options;        % Options for the solver
    hotStartQp;     % Store the QP for qpOASES hot start
    isColdStart;    % Boolean for qpOASES hot/cold start
end

methods
    function setup(obj)
    % Setup the solver
        obj.options = qpOASES_options('reliable');
        obj.isColdStart = true;
    end
    
    function [x,cost,exitflag] = solve(obj,H,f,A,b,lb,ub)
    % Solve a QP problem defined as follow
    % min 1/2*x'*H*x + f'*x   subject to:  A*x <= b  and lb <= x <= ub
        if obj.isColdStart
            [obj.hotStartQp,x,cost,exitflag] = qpOASES_sequence(...
                'i',H,f,A,lb,ub,[],b,obj.options...
            );
            obj.isColdStart= false;
        else
            try
                [x,cost,exitflag] = qpOASES_sequence(...
                    'm',obj.hotStartQp,H,f,A,lb,ub,[],b,obj.options...
                );
            catch ME
                if (strcmp(ME.message,'ERROR (qpOASES): Hotstart failed.'))
                    % Hot start failed: Delete previous QP and start a new
                    % one
                    warning('qpOASES: Hotstart failed.')
                    qpOASES_sequence('c',obj.hotStartQp);
                    [obj.hotStartQp,x,cost,exitflag] = qpOASES_sequence(...
                        'i',H,f,A,lb,ub,[],b,obj.options...
                    );
                else
                    rethrow(ME);
                end
            end
        end
    end
    
    function reset(~)
    % Reset the solver
        
    end
    
    function release(obj)
    % Release the solver
        if ~obj.isColdStart
            qpOASES_sequence('c',obj.hotStartQp);
        end
        obj.isColdStart = true;
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

