classdef DriverNMPC < matlab.System & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime
%DriverNMPC provide a MATLAB System to implement a Non-linear Model
%Predictive Controller as a driver model.
% 
% This class requires MPCTools and CasADi, make sure these packages are
% installed and added to the MATLAB path.
% 
% The driver aims to follow the path described by the road object. The MPC
% extracts waypoints from the road to create a spline that approximate the
% road path. This spline is then used in the cost function of the MPC to
% penalize the lateral position error and the heading error.
% 
% The spline is approximated using two polynomials, x(s) for the
% x-coordinates and y(s) for the y-coordinates. This also requires to
% compute the curvilinear abscissa s and therefore the road must be
% expressed using both a caartesian frame and a Frenet-Serret frame. This
% however helps to define maneuver when a simple polynomial y(x) could not
% be used (e.g. a U-turn).
% 
% The MPC uses an augmented plant for reference tracking, where the 
% augmented state is
%       z_k^T = [x_k^T u_k-1^T]
% and the new control input Delta u_k = u_k - u_k-1.
% 
% @author Louis Filipozzi

    % Public properties
    properties(Nontunable,Logical)
        % Show optimal cost
        showCost = false;
        % Show solving time
        showSolveTime = false;
        % Show optimal control sequence
        showControlSequence = false;
        % Show solver status
        showSolverStatus = false;
    end
    
    properties(Nontunable, PositiveInteger)
        % Number of preview point
        Nt = 10;
        % Spline order
        NSplineOrder = 4;
    end

    properties(Hidden, Constant, PositiveInteger)
        % Number of states
        Nx = 5;
        % Number of inputs
        Nu = 2;
        % Number of collocation point
        Nc = 2;
    end

    properties(Nontunable)
        % Preview time
        Tpreview = 2;
        % Road arc length
        sRoad = [0 1];
        % Road x coord.
        xRoad = [0 1];
        % Road y coord.
        yRoad = [0 0];
        % Wheelbase
        L = 2.78;
    end
    
    % State
    properties(DiscreteState)
        % Control inputs at the last timestep
        uprev
    end

    % Pre-computed constants
    properties(Nontunable, Access = private)
        solver  % MPCTools solver
        NxAug   % Number of state of the augmented model
        NuAug   % Number of input of the augmented model
        Ts      % Sampling time
    end

    methods(Access = protected)
        %% setupImpl & stepImpl & resetImpl
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            
            % Compute sampling time
            obj.Ts = obj.Tpreview / obj.Nt;
            
            % Compute sizes of the augmented model
            obj.NxAug = obj.Nx + obj.Nu;
            obj.NuAug = obj.Nu;
            
            % Build the solvers
            obj.solver = obj.buildMPCSolver();
        end

        function varargout = stepImpl(obj, varargin)
            % Implement algorithm. 
            
            % Setup chronometer
            if obj.showSolveTime
                tic;
            end
            
            % Solve MPC problem
            [du_seq, cost] = obj.solveMPC(varargin{:});
            
            % Compute the new control input
            du = du_seq(:,1);
            u = obj.uprev + du;
            
            % Save the new control input for the next timestep
            obj.uprev = u;
            
            % Return only the control input without the known inputs
            u_output = u(1:obj.Nu);
            
            % Returns output: [u_output, solveTime, cost]
            varargout{1} = u_output;
            if obj.showSolveTime
                solveTime = toc;
                varargout = [varargout{:}, {solveTime}];
            end
            if obj.showCost
                varargout = [varargout{:}, {cost}];
            end
            if obj.showControlSequence
                sequence = zeros(obj.Nu, obj.Nt);
                sequence(1:obj.Nu,1) = u_output;
                for k = 2:obj.Nt
                    sequence(1:obj.Nu,k) = ...
                        sequence(1:obj.Nu,k-1) + du_seq(1:obj.Nu,k);
                end
                varargout = [varargout{:}, {sequence}];
            end
            if obj.showSolverStatus
                varargout = [varargout{:}, ...
                    {double(isequal(obj.solver.status, ...
                    'Solve_Succeeded'))}];
            end
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.uprev = zeros(obj.NuAug,1);
        end
    end
    
    methods(Static, Access=private)
        %% Vehicle and road equation of motion
        function state_dot = MpcStateEquation(state, input, ...
                coefX, coefY, L)
            % Geometric model (with Frennet-Serret coordinates)
            % Unpack states
            U     = state(1);
            wz    = state(2);
            s     = state(3);
            sn    = state(4);
            alpha = state(5);
            
            % Unpack inputs
            delta = input(1);
            ax    = input(2);
            
            % Compute the road curvature using the x and y spline
            xd = 0; yd = 0; xdd = 0; ydd = 0;
            for k = 2:numel(coefX)
                xd = xd + (k-1)*coefX(k) * s.^(k-2);
            end
            for k = 2:numel(coefY)
                yd = yd + (k-1)*coefY(k) * s.^(k-2);
            end
            for k = 3:numel(coefX)
                xdd = xdd + (k-1)*(k-2)*coefX(k) * s.^(k-3);
            end
            for k = 3:numel(coefY)
                ydd = ydd + (k-1)*(k-2)*coefY(k) * s.^(k-3);
            end
            kappa = (xd*ydd - xdd*yd) ./ (xd.^2 + yd.^2).^(3/2);
            
            % Compute state derivatives (vehicle model)
            U_dot     = ax;
            wz_dot    = U/L*delta;
            psi_dot   = wz;
            s_dot     = U*cos(alpha)/(1-sn*kappa);
            sn_dot    = U*sin(alpha);
            alpha_dot = psi_dot - kappa*s_dot;
            
            % Stack states
            state_dot = [U_dot; wz_dot;
                s_dot; sn_dot; alpha_dot];
        end
        
        %% Stage cost
        function cost = MpcStageCost(z, du, Uref)
            % Define the stage cost (cost for one timestep of the
            % prediction horizon)
            
            trackError = z(4); % Cross track error (lateral position error)
            yawError   = z(5); % Heading angle error
            UError     = z(1) - Uref;   % Speed error
            e = [trackError; yawError; UError];
            
            % Weights
            qTrackError    =  50;
            qYawError      =  50 / (pi/4)^2;
            qUError        =  10;
            qdSteering     =   1 / (pi/4)^2;
            qdAcceleration =   1 / 9.81^2;
            
            % Weight matrices
            Qe  = diag([qTrackError, qYawError, qUError]);
            Rdu = diag([qdSteering, qdAcceleration]);
            
            cost = e'*Qe*e + du'*Rdu*du;
        end
        
        %% Terminal cost
        function cost = MpcTerminalCost(~)
            cost = 0;
        end
        
        %% Optimal problem constraints
        function ineq = MpcConstraints(z,du)
            ineq = [];
            % Constraint the steering angle
            delta = z(6) + du(1);
            deltaMax = pi/2;
            ineq = [ineq;
                delta - deltaMax;
               -delta - deltaMax];
            % Constraint the longitudinal acceleration
            acc = z(7) + du(2);
            accMax = 9.81;
            ineq = [ineq;
                acc - accMax;
               -acc - accMax];
        end
    end
    
    methods(Access=private)        
        %% Approximate road waypoints by a spline
        function [coefX, coefY] = getRoadSplineCoef(obj, U, sStart)
            % Reduce the number or points used to compute the splines
            NWaypoints = 2*obj.Nt;    % Number of waypoints
            
            % Make sure U is non-zero
            U = max(1,abs(U)) * sign(U);    
            
            % Compute the x and y coordinates of the waypoints
            sWaypoints = sStart + linspace(0,2*U*obj.Tpreview,NWaypoints);
            xWaypoints = interp1(obj.sRoad,obj.xRoad,sWaypoints',...
                'linear','extrap');
            yWaypoints = interp1(obj.sRoad,obj.yRoad,sWaypoints',...
                'linear','extrap');
            % Reset initial arc length to zero
            sWaypoints = sWaypoints - sStart;
            
            % Compute the spline coefficients
            sSample = ones(NWaypoints, obj.NSplineOrder+1);
            for k = 1:obj.NSplineOrder
                sSample(:,k+1) = sWaypoints.^k;
            end
            coefX = pinv(sSample)*xWaypoints;
            coefY = pinv(sSample)*yWaypoints;
        end
        
        %% Build MPC solver
        function solver = buildMPCSolver(obj)
            % This function is used to build the MPC solver.
            % This MPC has an augmented formulation in order to follow a
            % requested target.
            
            % Import MPC Tools
            mpcBuilder = import_mpctools();
            
            % Create structure to settings of the solver
            kwargsSolver = struct;
            
            % +---------------------------+
            % |   Define model and cost   |
            % +---------------------------+
            % Augment the model
            %         [x_k+1]   [ EOM(x_k, u_k-1 + Delta u-k, param) ]
            % z_k+1 = [ u_k ] = [          u_k-1 + Delta u-k         ]
            % Augmented state z_k^T = [x_k^T u_k-1^T]
            
            % Add custom variables
            customvar = {'coefX','coefY','Uref'};
            kwargsSolver.customvar = customvar;
            
            % Define state equation
            stateEquation = @(z,du,coefX,coefY) ...
                [obj.MpcStateEquation(...
                    z(1:obj.Nx), z((obj.Nx+1):(obj.Nx+obj.Nu))+du, ...
                    coefX, coefY, obj.L);
                z((obj.Nx+1):(obj.Nx+obj.Nu)) + du];
            fmpc = mpcBuilder.getCasadiFunc(stateEquation, ...
                [obj.NxAug, obj.NuAug, ...
                obj.NSplineOrder+1, obj.NSplineOrder+1], ...
                {'x', 'u','coefX','coefY'}, {'dint'});

            % Define stage cost
            stagecost = @(z,du,Uref) ...
                obj.MpcStageCost(z,du,Uref);
            l  = mpcBuilder.getCasadiFunc(stagecost, ...
                [obj.NxAug, obj.NuAug, 1], ...
                {'x','u','Uref'}, {'l'});

            % Define terminal cost
            termcost = @(z) obj.MpcTerminalCost(z);
            Vf = mpcBuilder.getCasadiFunc(termcost, ...
                [obj.NxAug], {'x'}, {'Vf'});

            % Define number of state, input and horizon
            N = struct('x', obj.NxAug, 'u', obj.NuAug, 't', obj.Nt,...
                'c', obj.Nc);

            % The MPC sampling time must be provided for non-linear MPC
            kwargsSolver.Delta = obj.Ts;
            
            % +---------------------------+
            % |  Define slew rate bounds  |
            % +---------------------------+
            % Note: The model has been augmented, therefore the bounds on u
            % correspond to bounds on Delta u_k and the bounds on x
            % correspond to bounds on z_k
            lb = struct();
            ub = struct();
            lb.u = -inf(obj.NuAug, obj.Nt);
            ub.u =  inf(obj.NuAug, obj.Nt);
            lb.x = -inf(obj.NxAug, obj.Nt+1);
            ub.x =  inf(obj.NxAug, obj.Nt+1);

            % +---------------------------+
            % |  Define ineq constraints  |
            % +---------------------------+
            % Inequality constraints
            constraintFunc = @(z,du) obj.MpcConstraints(z,du);
            eineq = mpcBuilder.getCasadiFunc(constraintFunc, ...
                [obj.NxAug, obj.NuAug], {'x', 'u'}, {'eineq'});
            
            % +---------------------------+
            % |        Build solver       |
            % +---------------------------+
            % Set no verbose mode
            kwargsSolver.verbosity = 0;
            
            % Build MPC solver
            solver = mpcBuilder.nmpc('f', fmpc, 'N', N, ...
                'l', l, 'Vf', Vf, 'lb', lb, 'ub', ub, 'e', eineq, ...
                '**', kwargsSolver,'solver','ipopt');
        end
        
        %% Solve MPC problem
        function [du_seq, cost] = solveMPC(obj, varargin)
            % Solve the MPC Problem for a given solver
            
            % Rename inputs
            ref = varargin{1};
            x   = varargin{2};
            
            % Compute the extended state for augmented formulation
            % z_k^T = [x_k^T u_k-1^T]
            xaug = [x; obj.uprev];
            
            % Set initial state to current state
            obj.solver.fixvar('x', 1, xaug);
            
            % Compute spline coefficient and give it to the solver
            U = x(1);
            sStart = x(3);
            [coefX, coefY] = getRoadSplineCoef(obj, U, sStart);
            obj.solver.fixvar('coefX', 1, coefX);
            obj.solver.fixvar('coefY', 1, coefY);
            obj.solver.fixvar('Uref',  1, ref(1));
            
            % Solve the MPC problem
            obj.solver.solve()

            % Return error if the problem has not been solved
            if ~isequal(obj.solver.status, 'Solve_Succeeded')
                msg = ['Solver failed! Status is ', obj.solver.status];
                warning(msg);
            end

            % Return the first control input
            du_seq = obj.solver.var.u;
            
            % Save the cost
            cost = obj.solver.obj;
        end
    end
    
    methods(Access=protected)
        %% Specify sampling time
        function sts = getSampleTimeImpl(obj)
            % Define sampling time
            sts = createSampleTime(obj,'Type','Discrete',...
              'SampleTime',obj.Tpreview / obj.Nt,'OffsetTime',0);
        end
        
        %% Number of input and output ports
        % Indicate number of input ports of the system object
        function num = getNumInputsImpl(~)
            % Two inputs minimum for state feedback and reference signal
            num = 2;
        end
        
        % Indicate the number of output ports of the system object
        function num = getNumOutputsImpl(obj)
            % One output minimum for the control input
            num = 1;
            if obj.showCost
                num = num + 1;
            end
            if obj.showSolveTime
                num = num + 1;
            end
            if obj.showControlSequence
                num = num + 1;
            end
            if obj.showSolverStatus
                num = num + 1;
            end
        end
        
        %% Inputs and outputs port name
        function varargout = getInputNamesImpl(~)
            % Non-optional inputs
            varargout{1} = 'reference';
            varargout{2} = 'state';
        end
   
        function varargout = getOutputNamesImpl(obj)
            % Non-optional outputs
            varargout{1} = 'control';
            % Optional outputs
            if obj.showSolveTime
                varargout = [varargout{:}, {'time'}];
            end
            if obj.showCost
                varargout = [varargout{:}, {'cost'}];
            end
            if obj.showControlSequence
                varargout = [varargout{:}, {'sequence'}];
            end
            if obj.showSolverStatus
                varargout = [varargout{:}, {'status'}];
            end
       end
        
        %% Input and output ports properties (size, data type, ...)
        function ds = getDiscreteStateImpl(~)
            % Return structure of properties with DiscreteState attribute
            ds = struct([]);
        end
            
        function flag = isInputSizeLockedImpl(~,~)
            % Inputs are not allowed to change in size
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
            
            if obj.showSolverStatus
                varargout = [varargout{:}, {true}];
            end
        end

        function varargout = getOutputSizeImpl(obj)
            
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
            
            % solver status
            if obj.showSolverStatus
                varargout = [varargout{:}, {[1 1]}];
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
            
            % solver status
            if obj.showSolverStatus
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
            
            % solver status
            if obj.showSolverStatus
                varargout = [varargout{:}, {false}];
            end
        end
        
        %% Customize MATLAB System block apperance
        function icon = getIconImpl(~)
            icon = {'NMPC','Driver Model'};
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
            header = matlab.system.display.Header('DriverNMPC',...
                'Title','NMPC Driver model',...
                'Text', ['NMPC implementation of a driver model for ', ...
                'longitudinal and lateral control with preview.']);
        end
        
        function groups = getPropertyGroupsImpl
            % Customize dialog box
            
            
            % +-----------------------+
            % |    MPC parameters     |
            % +-----------------------+
            % Create section for preview parameters
            previewSection = matlab.system.display.Section(...
                'Title','Preview', ...
                'PropertyList',{'Tpreview','Nt'});
            
            % Create section for road information
            roadSection = matlab.system.display.Section(...
                'Title','Road', ...
                'PropertyList',{'sRoad','xRoad','yRoad','NSplineOrder'});
            
            % Create section for vehicle information
            vehicleSection = matlab.system.display.Section(...
                'Title','Vehicle', ...
                'PropertyList',{'L'});
            
            % Create tab for MPC parameters configuration
            MPCParametersGroup = matlab.system.display.SectionGroup(...
                'Title','Parameters',...
                'Sections',[previewSection, roadSection, vehicleSection]);
            
            % +-----------------------+
            % |   MPC optional ports  |
            % +-----------------------+
            % Create section for optional output ports
            additionalOutportsSection = matlab.system.display.Section(...
                'Title','Additional Outports', ...
                'PropertyList',{'showCost','showSolveTime',...
                'showControlSequence', 'showSolverStatus'});
            
            % Create tab for general configuration
            generalGroup = matlab.system.display.SectionGroup(...
                'Title','Ports',...
                'Sections', additionalOutportsSection);
            
            % +-----------------------+
            % |      Create tabs      |
            % +-----------------------+
            groups = [MPCParametersGroup, generalGroup];
       end
    end
end
