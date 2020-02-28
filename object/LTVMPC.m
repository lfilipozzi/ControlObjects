classdef LTVMPC < matlab.System & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime
%LTVMPC provides an implementation of an MPC in Simulink. 
%   The MPC problem is solved using the batch approach with MATLAB's 
%   quadprog.
%
%   This MATLAB System Object support the following features:
%    - Known inputs in the MPC model
%    - Soft constraints
%    - Different control and prediction horizon
% 
% @author Louis Filipozzi

    properties(Nontunable,Logical)
        % Show optimal cost
        showCost = false;
        % Show solving time
        showSolveTime = false;
        % Show optimal control sequence
        showControlSequence = false;
        % Show slack variable
        showSlackVariable = false;
        % Known model inputs
        hasKnownInput = false;
        % Inequality constrainsts
        hasIneqConstraints = false;
        % Actuator bounds
        hasActuatorBounds = false;
        % Soft inequality constraints
        hasSoftConstraints = false;
        % Add scaling factors
        hasScaling = false;
    end
    
    properties(Nontunable, PositiveInteger)
        % Control horizon
        Nt = 10;
        % Prediction horizon
        Np = 10;
        % Number of states
        Nx = 7;
        % Number of inputs
        Nu = 4;
    end

    properties(Nontunable)
        % Number of known inputs
        Nw = 0;
        % Sampling time
        Ts = 0.01;
        % Indices
        softConstraintsIndex = {0};
        % Weight
        softConstraintsWeight = 1e6;
        % Actuator and known inputs
        actuatorScaling = 1;
        % Slack variables
        slackScaling = 1;
    end

    % Pre-computed constants
    properties(Access = private)
        % Options for the quadprog solver
        quadprog_option
        % Number of slack variable for soft constraints
        Nslack
    end

    methods(Access = protected)
        %% setupImpl & stepImpl & resetImpl
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            
            % Set quadprog options
            obj.quadprog_option = optimoptions('quadprog',...
                'display','off', ...
                'MaxIterations',1000, ...           % Default: 200
                'OptimalityTolerance', 1e-10, ...   % Default: 1e-8
                'StepTolerance', 1e-16, ...         % Default: 1e-12
                'ConstraintTolerance', 1e-8 ...     % Default: 1e-8
                );
            
            % If no known input, set Nw to 0 automatically
            if ~obj.hasKnownInput
                obj.Nw = 0;
            end
            
            % Check number of slack variables
            if obj.hasSoftConstraints
                obj.Nslack = numel(obj.softConstraintsIndex);
            else
                obj.Nslack = 0;
            end
        end

        function varargout = stepImpl(obj,varargin)
            % Implement algorithm.
            
            % Rename inputs
            Q  = varargin{1};
            R  = varargin{2};
            T  = varargin{3};
            fx = varargin{4};
            fu = varargin{5};
            x0 = varargin{6};
            A  = varargin{7};
            B  = varargin{8};
            
            num = 8;
            
            % Known inputs
            if obj.hasKnownInput
                num = num + 1;
                wmeas = varargin{num};
            else
                wmeas = [];
            end
            
            % Actuator saturation bounds
            if obj.hasActuatorBounds
                num = num + 2;
                ulb = [varargin{num - 1}; wmeas];
                uub = [varargin{num};     wmeas];
            else
                ulb = -inf(obj.Nu + obj.Nw,1);
                uub =  inf(obj.Nu + obj.Nw,1);
            end
            
            % Inequality constraints
            if obj.hasIneqConstraints
                num = num + 3;
                Aineqx = varargin{num-2};
                Ainequ = varargin{num-1};
                bineq  = varargin{num};
            else
                Aineqx = zeros(0, obj.Nx);
                Ainequ = zeros(0, obj.Nu + obj.Nw);
                bineq  = [];
            end
            Nineq = size(bineq,1);
            Aineqs = zeros(Nineq, obj.Nslack);
            for i = 1:obj.Nslack
                Aineqs(obj.softConstraintsIndex{i},i) = -1;
            end
            
            % MPC problem formulation
            MPC_QP = struct('A',A,'B',B,'x0',x0,...
                'Q',Q,'R',R,'T',T,'fx',fx,'fu',fu,...
                'ulb',ulb,'uub',uub,...
                'Aineqx',Aineqx,'Ainequ',Ainequ,...
                'Aineqs',Aineqs,'bineq',bineq...
            );
            
            % Setup chronometer
            if obj.showSolveTime
                tic;
            end
            
            % Solve MPC problem
            [u_opt, cost, exitflag, slack] = obj.solveMPC(MPC_QP);
            u_opt = reshape(u_opt, [obj.Nu+obj.Nw obj.Nt]);
            
            % Warning if the problem has not been solved
            if exitflag < 1
                errMessage = obj.getWarningMessage(exitflag);
                warning(errMessage);
            end
            
            % Record time to solve
            if obj.showSolveTime
                solveTime = toc;
            end
            
            % Returns output: [u, solveTime, cost, exitflag]
            varargout{1} = u_opt(1:obj.Nu,1);
            if obj.showSolveTime
                varargout = [varargout{:}, {solveTime}];
            end
            if obj.showCost
                varargout = [varargout{:}, {cost}];
            end
            if obj.showControlSequence
                varargout = [varargout{:}, {u_opt(1:obj.Nu,:)}];
            end
            varargout = [varargout{:}, {exitflag}];
            if obj.hasSoftConstraints && obj.showSlackVariable
                varargout = [varargout{:}, {slack}];
            end
        end

        function resetImpl(~)
            % Initialize / reset discrete-state properties
            
        end
        
        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values
            
            % Check the prediction horizon is at least as long as the
            % control horizon
            if obj.Np < obj.Nt
                msg = ['The prediction horizon must be larger than the '...
                    'control horizon'];
                error(msg)
            end
        end
    end
    
    methods(Access = private)
        %% Solve MPC problem
        function [U, cost, exitflag, slack] = solveMPC(obj,MPC_QP)
            % Solve the LTV MPC problem using batch approach.
            
            % Unpack matrices
            Q   = MPC_QP.Q;
            R   = MPC_QP.R;
            T   = MPC_QP.T;
            fx  = MPC_QP.fx;
            fu  = MPC_QP.fu;
            A   = MPC_QP.A;
            B   = MPC_QP.B;
            x0  = MPC_QP.x0;
            ulb = MPC_QP.ulb;
            uub = MPC_QP.uub;
            Aineqx = MPC_QP.Aineqx;
            Ainequ = MPC_QP.Ainequ;
            Aineqs = MPC_QP.Aineqs;
            bineq   = MPC_QP.bineq;
            
            % Get the number of states, inputs, constraints and horizon
            Nx_MPC = size(A,1);         % Number of states
            Nu_MPC = size(B,2);         % Number of control inputs
            Nt_MPC = obj.Nt;            % Control horizon
            Np_MPC = obj.Np;            % Prediction horizon
            Nineq  = size(Aineqx, 1);   % Ineqaulity constraints
            
            % +-------------------------+
            % |      Cost function      |
            % +-------------------------+
            
            % Solve the problem by using a batch approch:
            % Write the state equation over the prediction horizon:
            %       x = Sx x0 + Su u
            % where:
            %       x = [x0 ... xN]^T
            %       u = [u0 ... uN-1]^T
            
            Sx = zeros(Nx_MPC*Np_MPC, Nx_MPC);
            Sx(1:Nx_MPC, 1:Nx_MPC) = eye(Nx_MPC);
            for i = 1:Np_MPC-1
                Sx(i*Nx_MPC + (1:Nx_MPC), 1:Nx_MPC) = ...
                    A * Sx((i-1)*Nx_MPC + (1:Nx_MPC), 1:Nx_MPC);
            end

            Su = zeros(Nx_MPC*Np_MPC,Nu_MPC*Nt_MPC);
            Su(Nx_MPC+(1:Nx_MPC), 1:Nu_MPC) = B;
            % Fill first column
            for i = 2:Np_MPC-1
                Su(Nx_MPC*i+(1:Nx_MPC), 1:Nu_MPC) = ...
                    A * Su(Nx_MPC*(i-1)+(1:Nx_MPC), 1:Nu_MPC);
            end
            % Copy first column to other columns
            for j = 1:Nt_MPC-1
                Su(Nx_MPC*(j+1)+1:end, Nu_MPC*j+(1:Nu_MPC)) = ...
                    Su(Nx_MPC+1:Nx_MPC*(Np_MPC-j), 1:Nu_MPC);
            end
            
            % Write cost function for batch approach
            QQ = zeros(Nx_MPC*Np_MPC);
            for k = 1:Np_MPC
                QQ((k-1)*Nx_MPC+1:k*Nx_MPC,(k-1)*Nx_MPC+1:k*Nx_MPC) = Q;
            end

            RR = zeros(Nu_MPC*Nt_MPC);
            for k = 1:Nt_MPC
                RR((k-1)*Nu_MPC+1:k*Nu_MPC,(k-1)*Nu_MPC+1:k*Nu_MPC) = R;
            end
            
            TT = zeros(Nx_MPC*Np_MPC,Nu_MPC*Nt_MPC);
            for k = 1:Nt_MPC
                TT((k-1)*Nx_MPC+1:k*Nx_MPC,(k-1)*Nu_MPC+1:k*Nu_MPC) = T;
            end
            
            Fx = zeros(Nx_MPC*Np_MPC, 1);
            for k = 1:Np_MPC
                Fx((k-1)*Nx_MPC+1:k*Nx_MPC, 1) = fx;
            end
            
            Fu = zeros(Nu_MPC*Nt_MPC, 1);
            for k = 1:Nt_MPC
                Fu((k-1)*Nu_MPC+1:k*Nu_MPC, 1) = fu;
            end
            
            H = 2 * (Su' * QQ * Su + TT' * Su + Su' * TT + RR);
            f = Su'*Fx + Fu + 2 * (QQ * Su + TT)' * Sx * x0;
            
            % +-------------------------+
            % |         Bounds          |
            % +-------------------------+
            % Check min <= max
            if sum(ulb > uub)
                error('Unfeasible bounds on u')
            end
            
            % Extend matrices dumin, dumax to obtain dulb and duub where:
            %       Ulb <= U <= Uub
            % with:
            %       U = [u0; u1; ... uN-1]
            Ulb = zeros(Nu_MPC*Nt_MPC,1);
            for k = 1:Nt_MPC
                Ulb((k-1)*Nu_MPC+1:k*Nu_MPC,:) = ulb;
            end
            
            Uub = zeros(Nu_MPC*Nt_MPC,1);
            for k = 1:Nt_MPC
                Uub((k-1)*Nu_MPC+1:k*Nu_MPC,:) = uub;
            end
            
            % +-------------------------+
            % | Inequality  constraints |
            % +-------------------------+
            % Inequalities are written as:
            %       Aineqx x + Ainequ u <= bineq
            
            if isempty(Aineqx) && isempty(Ainequ)
                AAineq = zeros(0,Nu_MPC*Nt_MPC);
                Bineq = [];
            else
                % Create matrix blockdiag(Au,...,Au)
                AAinequ = zeros(Nineq*Np_MPC, Nu_MPC*Nt_MPC);
                for k = 1:Nt_MPC
                    AAinequ((k-1)*Nineq+1:k*Nineq,...
                        (k-1)*Nu_MPC+1:k*Nu_MPC) = Ainequ;
                end
                
                % Create matrix blockdiag(Ax,...,Ax)
                AAineqx = zeros(Nineq*Np_MPC, Nx_MPC*Np_MPC);
                for k = 1:Np_MPC
                    AAineqx((k-1)*Nineq+1:k*Nineq,...
                        (k-1)*Nx_MPC+1:k*Nx_MPC) = Aineqx;
                end
                
                % Write inequalities as: AAineq * U <= Bineq
                AAineq = AAineqx * Su + AAinequ;
                Bineq = zeros(Nineq*Np_MPC,1);
                for k = 1:Np_MPC
                    Bineq((k-1)*Nineq+1:k*Nineq) = ...
                        bineq;
                end
                Bineq = Bineq - AAineqx * Sx * x0;
            end
            
            % +-----------------------------------+
            % |    Soft inequality constraints    |
            % +-----------------------------------+
            if obj.hasSoftConstraints
                % Compute the augmented matrix 
                %   [Aineq_slack]
                %   [     ...   ]
                %   [Aineq_slack]
                % for inequality constraints
                AAineqs = zeros(size(AAineq,1), obj.Nslack);
                for k = 1:Np_MPC
                    AAineqs((k-1)*size(Aineqs,1)+1:...
                        k*size(Aineqs,1),:) = ...
                        Aineqs;
                end
                
                % Modify the optimization problem
                H = blkdiag(H, diag(obj.softConstraintsWeight));
                f = [f; zeros(obj.Nslack, size(f,2))];
                AAineq = [AAineq AAineqs];
                Ulb = [Ulb;  zeros(obj.Nslack,1)];
                Uub = [Uub;  inf(obj.Nslack,1)];
            end
            
            % +-------------------------+
            % |    Scale the problem    |
            % +-------------------------+
            % Scale the QP problem
            if obj.hasScaling
                Ddiag = ones(size(H,1),1);
                for k = 1:Nt_MPC
                    Ddiag((k-1)*Nu_MPC+(1:6)) = obj.actuatorScaling;
                end
                if obj.hasSoftConstraints
                    Ddiag(end-obj.Nslack+1:end) = obj.slackScaling;
                end
                Ddiag = 1./Ddiag;
                D = diag(Ddiag);
                H = D'*H*D;
                f = Ddiag .* f;
                AAineq = AAineq*D;
                Ulb  = Ddiag .\ Ulb;
                Uub  = Ddiag .\ Uub;
            end
            
            % +-------------------------+
            % |      Solve problem      |
            % +-------------------------+
            H = (H+H')/2;   % Make sure H is symmetric
            [U,~,exitflag] = quadprog(H,f,...
                AAineq,Bineq,[],[],Ulb,Uub,...
                [],... 
                obj.quadprog_option);
            
            if isempty(U)
                error('Solver failed! No solution.')
            end
            
            % +-------------------------+
            % |      Remove scaling     |
            % +-------------------------+
            if obj.hasScaling
                U = Ddiag .* U;
            end
            
            % +-------------------------+
            % |      Compute output     |
            % +-------------------------+
            slack = [];
            if obj.hasSoftConstraints
                % Separate the slack variables from the control inputs
                slack = U(end-obj.Nslack+1:end);
                U = U(1:end-obj.Nslack);
            end
            
            % Compute the state over the prediction horizon
            X = Sx * x0 + Su * U;
            
            % Compute the cost (no constant term, no soft constraints)
            cost = X' * QQ * X + U' * RR * U + ...
                X' * TT * U + U' * TT' * X;
        end
        
        %% quadprog's custom error message
        function errMessage = getWarningMessage(obj,exitflag)
            errMessage = 'Solver failed at time ';
            currentTimestep = obj.getCurrentTime();
            errMessage = [errMessage, num2str(currentTimestep,10), '! '];
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
    
    methods(Access=protected)
        %% Specify sampling time
        function sts = getSampleTimeImpl(obj)
            % Define sampling time
            sts = createSampleTime(obj,'Type','Discrete',...
              'SampleTime',obj.Ts,'OffsetTime',0);
        end
        
        %% Number of input and output ports
        % Indicate number of input ports of the system object
        function num = getNumInputsImpl(obj)
            % Two inputs minimum for state feedback and reference signal
            num = 8;
            if obj.hasKnownInput
                % Add inport for known input
                num = num + 1;
            end
            if obj.hasActuatorBounds
                % Add two signals for max and min actuator bounds
                num = num + 2;
            end
            if obj.hasIneqConstraints
                % Add one input for ineqaulity bounds
                num = num + 3;
            end
        end
        
        % Indicate the number of output ports of the system object
        function num = getNumOutputsImpl(obj)
            % Two outputs minimum for the control input and QP exitflag
            num = 2;
            if obj.showCost
                num = num + 1;
            end
            if obj.showSolveTime
                num = num + 1;
            end
            if obj.showControlSequence
                num = num + 1;
            end
            if obj.hasSoftConstraints && obj.showSlackVariable
                num = num + 1;
            end
        end
        
        %% Inputs and outputs port name
        function varargout = getInputNamesImpl(obj)
            % Non-optional inputs
            varargout{1} = 'Q';
            varargout{2} = 'R';
            varargout{3} = 'T';
            varargout{4} = 'fx';
            varargout{5} = 'fu';
            varargout{6} = 'x0';
            varargout{7} = 'A';
            varargout{8} = 'B';
            % Optional inputs
            if obj.hasKnownInput
                varargout = [varargout{:}, {'uknwown'}];
            end
            if obj.hasActuatorBounds
                varargout = [varargout{:}, {'ulb','uub'}];
            end
            if obj.hasIneqConstraints
                varargout = [varargout{:}, {'Aineqx','Ainequ','bineq'}];
            end
            if obj.hasIneqConstraints && obj.hasSoftConstraints
                varargout = [varargout{:}, {'Aineqs'}];
            end
        end
   
        function varargout = getOutputNamesImpl(obj)
            % Non-optional outputs
            varargout{1} = 'u_opt';
            % Optional outputs
            if obj.showSolveTime
                varargout = [varargout{:}, {'time'}];
            end
            if obj.showCost
                varargout = [varargout{:}, {'cost'}];
            end
            if obj.showControlSequence
                varargout = [varargout{:}, {'U_opt'}];
            end
            varargout = [varargout{:}, {'status'}];
            if obj.hasSoftConstraints && obj.showSlackVariable
                varargout = [varargout{:}, {'slack'}];
            end
       end
        
        %% Simulink functions
        function ds = getDiscreteStateImpl(~)
            % Return structure of properties with DiscreteState attribute
            ds = struct([]);
        end

        function flag = isInputSizeLockedImpl(~,~)
            % Input size is not allowed to change
            flag = true;
        end
        
        function varargout = isOutputFixedSizeImpl(obj)
            % Define if the output size if fixed
            
            % u_output
            varargout{1} = true;
            
            % solving time
            if obj.showSolveTime
                varargout = [varargout{:}, {true}];
            end
            
            % optimal cost
            if obj.showCost
                varargout = [varargout{:}, {true}];
            end
            
            % control sequence
            if obj.showControlSequence
                varargout = [varargout{:}, {true}];
            end
            
            % exitflag
            varargout = [varargout{:}, {true}];
            
            % slack variable
            if obj.hasSoftConstraints && obj.showSlackVariable
                varargout = [varargout{:}, {true}];
            end
        end

        function varargout = getOutputSizeImpl(obj)
            % Return size for each output port
            
            varargout{1} = [obj.Nu 1];
            
            % solving time
            if obj.showSolveTime
                varargout = [varargout{:}, {[1 1]}];
            end
            
            % optimal cost
            if obj.showCost
                varargout = [varargout{:}, {[1 1]}];
            end
            
            % control sequence
            if obj.showControlSequence
                varargout = [varargout{:}, {[obj.Nu obj.Nt]}];
            end
            
            % exitflag
            varargout = [varargout{:}, {[1 1]}];
            
            % slack variables
            if obj.hasSoftConstraints && obj.showSlackVariable
                varargout = [varargout{:}, ...
                    {[numel(obj.softConstraintsIndex) 1]}];
            end
        end
        
        function varargout = getOutputDataTypeImpl(obj)
            % u_output
            varargout{1} = 'double';
            
            % solving time
            if obj.showSolveTime
                varargout = [varargout{:}, {'double'}];
            end
            
            % optimal cost
            if obj.showCost
                varargout = [varargout{:}, {'double'}];
            end
            
            % control sequence
            if obj.showControlSequence
                varargout = [varargout{:}, {'double'}];
            end
            
            % exitflag
            varargout = [varargout{:}, {'double'}];
            
            % slack variables
            if obj.hasSoftConstraints && obj.showSlackVariable
                varargout = [varargout{:}, {'double'}];
            end
        end
        
        function varargout = isOutputComplexImpl(obj)
            % u_output
            varargout{1} = false;
            
            % solving time
            if obj.showSolveTime
                varargout = [varargout{:}, {false}];
            end
            
            % optimal cost
            if obj.showCost
                varargout = [varargout{:}, {false}];
            end
            
            % control sequence
            if obj.showControlSequence
                varargout = [varargout{:}, {false}];
            end
            
            % exitflag
            varargout = [varargout{:}, {false}];
            
            % slack variable
            if obj.hasSoftConstraints && obj.showSlackVariable
                varargout = [varargout{:}, {false}];
            end
        end
        
        %% Set inactive properties
        function flag = isInactivePropertyImpl(obj,propertyName)
            % Nw is needed only if we have known inputs
            if ismember(propertyName, {'Nw'})
                flag = not(obj.hasKnownInput);
            % Soft constraints properties only with soft constraints
            elseif ismember(propertyName, {'softConstraintsIndex',...
                    'softConstraintsWeight','showSlackVariable'})
                flag = not(obj.hasSoftConstraints);
            % Scaling used for QP solver
            elseif ismember(propertyName, {'actuatorScaling'})
                flag = not(obj.hasScaling);
            elseif ismember(propertyName, {'slackScaling'})
                flag = not(obj.hasScaling && obj.hasSoftConstraints);
            % All other properties are active
            else
                flag = false;
            end
        end
        
        %% Customize MATLAB System block apperance
        function icon = getIconImpl(~)
            icon = {'Linear Time','Varying','MPC'};
        end
    end
    
    %% Customization function of the MATLAB System
    methods (Static, Access = protected)
        function simMode = getSimulateUsingImpl
            % Force Interpreted Execution: Code generation is not supported
            simMode = 'Interpreted execution';
        end

        function flag = showSimulateUsingImpl
            % Return false if simulation mode hidden in System block dialog
            flag = false;
        end
        
        function isVisible = showFiSettingsImpl
            % Do not show Data type tab
            isVisible = false; 
        end
        
        function header = getHeaderImpl
            header = matlab.system.display.Header('Title','LTV MPC',...
                'Text','MATLAB System implementation of an LTV-MPC.');
        end
        
        function groups = getPropertyGroupsImpl
            % Customize dialog box
            
            % +-----------------------+
            % |    MPC Formulation    |
            % +-----------------------+
            % Create section for model
            MPCConstraints = matlab.system.display.Section(...
                'Title','Constraints',...
                'PropertyList',{'hasActuatorBounds',...
                'hasIneqConstraints',...
                'hasSoftConstraints','softConstraintsIndex',...
                'softConstraintsWeight'});
            
            % Create tab for MPC formulation
            MPCFormulationGroup = matlab.system.display.SectionGroup(...
                'Title','Formulation',...
                'Sections',[MPCConstraints]);
            
            % +-----------------------+
            % |    MPC parameters     |
            % +-----------------------+
            % Create section for hyper-parameters
            hyperParametersSection = matlab.system.display.Section(...
                'Title','Hyper-parameters', ...
                'PropertyList',{'Ts','Nt','Np'});
            
            % Create section for sizes
            sizesSection = matlab.system.display.Section(...
                'Title','Signals Sizes', ...
                'PropertyList',{'Nx','Nu'});
            
            % Create section for scaling the QP problem
            scalingSection = matlab.system.display.Section(...
                'Title','Scaling the Quadratic Problem', ...
                'PropertyList',{'hasScaling','actuatorScaling',...
                'slackScaling'});
            
            % Create tab for MPC parameters configuration
            MPCParametersGroup = matlab.system.display.SectionGroup(...
                'Title','Parameters',...
                'Sections',[hyperParametersSection, sizesSection, ...
                scalingSection]);
            
            % +-----------------------+
            % |   MPC optional ports  |
            % +-----------------------+
            % Create section for optional input ports
            additionalInportsSection = matlab.system.display.Section(...
                'Title','Additional Inports', ...
                'PropertyList',{'hasKnownInput','Nw'});
            
            % Create section for optional output ports
            additionalOutportsSection = matlab.system.display.Section(...
                'Title','Additional Outports', ...
                'PropertyList',{'showCost','showSolveTime',...
                'showControlSequence','showSlackVariable'});
            
            % Create tab for general configuration
            generalGroup = matlab.system.display.SectionGroup(...
                'Title','Ports',...
                'Sections', [additionalInportsSection,...
                additionalOutportsSection]);
            
            % +-----------------------+
            % |      Create tabs      |
            % +-----------------------+
            groups = [MPCFormulationGroup, MPCParametersGroup, ...
                generalGroup];
       end
        
    end
end
