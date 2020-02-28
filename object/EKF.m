classdef EKF < matlab.System & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime
% EKF Extended Kalman Filter
%
% Estimates the states of a discrete-time system.
%
% @author Louis Filipozzi
    
    % Public, nontunable, logical properties
    properties(Nontunable, Logical)
        % Add input port u
        hasInput = true;
        % State jacobian
        hasStateJacobianFcn = true;
        % Output jacobian
        hasOutputJacobianFcn = true;
        % Time-invariant Q
        isQTimeInvariant = true;
        % Time-invariant R
        isRTimeInvariant = true;
        % Output estimated model state x
        showEstimatedState = true;
        % Output state estimation covariance P
        showStateCovariance = false;
    end

    % Public, nontunable properties
    properties(Nontunable)
        % State transition function
        stateTransitionFcn = @(x,u) 0;
        % Jacobian
        stateJacobianFcn = @(x,u) 0;
        % Output function
        outputFcn = @(x,u) 0;
        % Jacobian
        outputJacobianFcn = @(x,u) 0;
        % Q
        QMat = 1;
        % R
        RMat = 1;
        % Initial state
        xInit = 1;
        % Initial state estimation error covariance
        PInit = 1;
        % Sampling time (-1 for inherited)
        Ts = 0.01;
    end
    
    % Public, tunable, positive integre properties
    properties(PositiveInteger)
        % Number of states
        NxUserDefined = 1;
        % Number of measurements
        NyUserDefined = 1;
        % Number of known inputs
        NuUserDefined = 1;
    end
    
    % Discrete states
    properties(DiscreteState)
        xHat;
        PHat;
        prevInput;
    end

    % Pre-computed positive-integer constants
    properties(Access = private, PositiveInteger)
        Nx;
        Ny;
        Nu;
    end

    methods(Access = protected)
        %% setupImpl & stepImpl & resetImpl
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            
            % Define size of the system
            [obj.Nx, obj.Nu, obj.Ny] = getSystemSize(obj);
        end

        function varargout = stepImpl(obj, varargin)
            % Rename input
            num = 0;
            % Add input for known inputs
            if obj.hasInput
                num = num + 1;
                uCurr = varargin{num};
                uPrev = obj.prevInput;
            else
                uCurr = 0;
                uPrev = 0;
            end
            % Add input for measurement
            num = num + 1;
            yCurr = varargin{num};

            % Add input for Q
            if obj.isQTimeInvariant
                Q = obj.QMat;
            else
                num = num + 1;
                Q = varargin{num};
            end
            
            % Add input for R
            if obj.isRTimeInvariant
                R = obj.RMat;
            else
                num = num + 1;
                R = varargin{num};
            end
            
            % Rename variables
            xPrevGPrev = obj.xHat;
            PPrevGPrev = obj.PHat;
            
            % A priori update
            [xCurrGPrev, PCurrGPrev] = obj.aPrioriUpdate(...
                Q, xPrevGPrev, PPrevGPrev, uPrev...
            );
            
            % A posteriori update
            [xCurrGCurr, PCurrGCurr] = obj.aPosterioriUpdate(...
                R, xCurrGPrev, PCurrGPrev, uCurr, yCurr...
            );
            
            % Store new states
            obj.xHat = xCurrGCurr;
            obj.PHat = PCurrGCurr;
            obj.prevInput = uCurr;
            
            % Return outputs
            varargout = {};
            % Output estimated model state x
            if obj.showEstimatedState
                varargout = [varargout{:}, {xCurrGCurr}];
            end
            % Output state estimation covariance P
            if obj.showStateCovariance
                varargout = [varargout{:}, {PCurrGCurr}];
            end
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.xHat = obj.xInit;
            obj.PHat = obj.PInit;
            obj.prevInput = zeros(obj.Nu,1);
        end

        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values
            
            % Check sampling time property
            if ~(obj.Ts == -1 || obj.Ts > 0)
                error('The sampling time must be positive or -1')
            end
        end
    end
    
    methods(Access = private)
        function [Nx, Nu, Ny] = getSystemSize(obj)
            Nx = obj.NxUserDefined;
            Ny = obj.NyUserDefined;
            if obj.hasInput
                Nu = obj.NuUserDefined;
            else
                Nu = 0;
            end
        end
        
        function [xCurrGPrev, PCurrGPrev] = aPrioriUpdate(...
                obj, Q, xPrevGPrev, PPrevGPrev, uPrev...
            )
        
            % Compute jacobian matrices
            A = obj.computeStateJacobian(xPrevGPrev, uPrev);
            
            % A priori update
            xCurrGPrev = obj.stateTransitionFcn(xPrevGPrev, uPrev);
            PCurrGPrev = A * PPrevGPrev * A' + Q;
        end
        
        function [xCurrGCurr, PCurrGCurr] = aPosterioriUpdate(...
                obj, R, xCurrGPrev, PCurrGPrev, uCurr, yCurr...
            )
            
            % Compute jacobian matrices
            C = obj.computeOutputJacobian(xCurrGPrev, uCurr);
            
            % Filter update
            Rtilde = C * PCurrGPrev * C' + R;
            L = PCurrGPrev * C' / Rtilde;
            
            % A posteriori update
            xCurrGCurr = xCurrGPrev + L * (...
                yCurr - obj.outputFcn(xCurrGPrev, uCurr) ...
            );
            PCurrGCurr = PCurrGPrev + - L * C * PCurrGPrev;
        end
        
        function A = computeStateJacobian(obj, xop, uop)
            % Compute the jacobian matrix of the state update equation.
            
            if obj.hasStateJacobianFcn
                % Compute jacobian matrices using user function
                A = obj.stateJacobianFcn(xop, uop);
            else
                % Compute jacobian matries using jaconbianest function
                nx = numel(xop);
                nu = numel(uop);
                updateFcn = @(x) obj.stateTransitionFcn(...
                    x(1:nx), x(nx+(1:nu))...
                );
                jacMat = jacobianest(updateFcn, [xop; uop]);
                A = jacMat(:,1:nx);
            end
        end
        
        function C = computeOutputJacobian(obj, xop, uop)
            % Compute the jacobian of the output equation.
            
            if obj.hasOutputJacobianFcn
                % Compute jacobian matrices using user function
                C = obj.outputJacobianFcn(xop, uop);
            else
                % Compute jacobian matries using jaconbianest function
                nx = numel(xop);
                nu = numel(uop);
                updateFcn = @(x) obj.outputFcn(...
                    x(1:nx), x(nx+(1:nu))...
                );
                jacMat = jacobianest(updateFcn, [xop; uop]);
                C = jacMat(:,1:nx);
            end
        end
    end
    
    methods(Access = protected)
        %% Specify sampling time
        function sts = getSampleTimeImpl(obj)
            % Define sampling time
            if obj.Ts == -1     % Inherited sampling time
                sts = createSampleTime(obj,'Type','Inherited');
            elseif obj.Ts > 0
                sts = createSampleTime(obj,'Type','Discrete',...
                  'SampleTime',obj.Ts,'OffsetTime',0);
            end
        end
        
        %% Number of input and output ports
        % Indicate number of input ports of the system object
        function num = getNumInputsImpl(obj)
            num = 0;
            % Add input for known inputs
            if obj.hasInput
                num = num + 1;
            end
            % Add input for measurement
            num = num + 1;
            % Add input for Q
            if ~obj.isQTimeInvariant
                num = num + 1;
            end
            % Add input for R
            if ~obj.isRTimeInvariant
                num = num + 1;
            end
        end
        
        % Indicate the number of output ports of the system object
        function num = getNumOutputsImpl(obj)
            num = 0;
            % Output estimated model state x
            if obj.showEstimatedState
                num = num + 1;
            end
            % Output state estimation covariance P
            if obj.showStateCovariance
                num = num + 1;
            end
        end
        
        %% Inputs and outputs port name
        function varargout = getInputNamesImpl(obj)
            varargout = {};
            % Add input for known inputs
            if obj.hasInput
                varargout = [varargout{:}, {'u'}];
            end
            % Add input for measurement
            varargout = [varargout{:}, {'y'}];
            % Add input for Q
            if ~obj.isQTimeInvariant
                varargout = [varargout{:}, {'Q'}];
            end
            % Add input for R
            if ~obj.isRTimeInvariant
                varargout = [varargout{:}, {'R'}];
            end
        end
   
        function varargout = getOutputNamesImpl(obj)
            varargout = {};
            % Output estimated model state x
            if obj.showEstimatedState
                varargout = [varargout{:}, {'xhat'}];
            end
            % Output state estimation covariance P
            if obj.showStateCovariance
                varargout = [varargout{:}, {'P'}];
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
            
            varargout = {};
            % Output estimated model state x
            if obj.showEstimatedState
                varargout = [varargout{:}, {true}];
            end
            % Output state estimation covariance P
            if obj.showStateCovariance
                varargout = [varargout{:}, {true}];
            end
        end

        function varargout = getOutputSizeImpl(obj)
            % Return size for each output port
            
            [Nbx, ~, ~] = getSystemSize(obj);
            
            varargout = {};
            % Output estimated model state x
            if obj.showEstimatedState
                varargout = [varargout{:}, {[Nbx, 1]}];
            end
            % Output state estimation covariance P
            if obj.showStateCovariance
                varargout = [varargout{:}, {[Nbx, Nbx]}];
            end
        end
        
        function varargout = getOutputDataTypeImpl(obj)
            varargout = {};
            % Output estimated model state x
            if obj.showEstimatedState
                varargout = [varargout{:}, {'double'}];
            end
            % Output state estimation covariance P
            if obj.showStateCovariance
                varargout = [varargout{:}, {'double'}];
            end
        end
        
        function varargout = isOutputComplexImpl(obj)
            varargout = {};
            % Output estimated model state x
            if obj.showEstimatedState
                varargout = [varargout{:}, {false}];
            end
            % Output state estimation covariance P
            if obj.showStateCovariance
                varargout = [varargout{:}, {false}];
            end
        end
        
        %% Set inactive properties
        function flag = isInactivePropertyImpl(obj,propertyName)
            
            % stateJacobianFcn required if jacobian is provided
            if ismember(propertyName, {'stateJacobianFcn'})
                flag = not(obj.hasStateJacobianFcn);
            % outputJacobianFcn required if jacobian is provided
            elseif ismember(propertyName, {'outputJacobianFcn'})
                flag = not(obj.hasOutputJacobianFcn);
            % QMat required only if time invariant
            elseif ismember(propertyName, {'QMat'})
                flag = not(obj.isQTimeInvariant);
            % RMat required only if time invariant
            elseif ismember(propertyName, {'RMat'})
                flag = not(obj.isRTimeInvariant);
            % Properties are active by default
            else
                flag = false;
            end
        end
        
        %% Customize MATLAB System block apperance
        function icon = getIconImpl(~)
            icon = {'Extended','Kalman','Filter'};
        end
    end
    
    %% Customization function of the MATLAB System
    methods (Static, Access = protected)
        function simMode = getSimulateUsingImpl
            simMode = "Interpreted execution";
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
            header = matlab.system.display.Header('EKF',...
                'Title','Extended Kalman Filter',...
                'Text', 'Estimates the states of a discrete-time system.');
        end

        function groups = getPropertyGroupsImpl
            % Customize dialog box
            
            % +-----------------------+
            % |    Model Parameters   |
            % +-----------------------+
            % Create section for system model
            SystemModelSection = matlab.system.display.Section(...
                'Title','System model',...
                'PropertyList',{...
                'stateTransitionFcn',...
                'hasStateJacobianFcn','stateJacobianFcn',...
                'outputFcn',...
                'hasOutputJacobianFcn','outputJacobianFcn',...
                'NxUserDefined','NuUserDefined','NyUserDefined',...
                'Ts'});
            
            % Create section for initial values
            InitialEstimateSection = matlab.system.display.Section(...
                'Title','Initial Estimates',...
                'PropertyList',{'xInit','PInit'});
            
            % Create section for noise characteristics
            NoiseCaracteristicsSection = matlab.system.display.Section(...
                'Title','Noise Characteristics',...
                'PropertyList',{'isQTimeInvariant','QMat',...
                'isRTimeInvariant','RMat'});
            
            % Create tab for MPC formulation
            ModelParametersGroup = matlab.system.display.SectionGroup(...
                'Title','Model Parameters',...
                'Sections',[SystemModelSection, InitialEstimateSection, ...
                NoiseCaracteristicsSection]);
            
            % +-----------------------+
            % |        Options        |
            % +-----------------------+
            % Create section for optional input ports
            additionalInportsSection = matlab.system.display.Section(...
                'Title','Additional Inports', ...
                'PropertyList',{'hasInput'});
            
            % Create section for optional output ports
            additionalOutportsSection = matlab.system.display.Section(...
                'Title','Additional Outports', ...
                'PropertyList',{'showEstimatedState',...
                'showStateCovariance'});
            
            % Create tab for general configuration
            OptionsGroup = matlab.system.display.SectionGroup(...
                'Title','Options',...
                'Sections', [additionalInportsSection,...
                additionalOutportsSection]);
            
            % +-----------------------+
            % |      Create tabs      |
            % +-----------------------+
            groups = [ModelParametersGroup, OptionsGroup];
       end
    end
end
