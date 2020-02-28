classdef UMV < matlab.System & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime
% UMV Unbiased Minimum Variance filter
%
% Estimates the states and unknown inputs of a strongly observable 
% discrete-time system.
%
% @author Louis Filipozzi
    
    % Public, nontunable, logical properties
    properties(Nontunable, Logical)
        % Add input port u
        hasKnownInput = true;
        % Time-invariant Q
        isQTimeInvariant = true;
        % Time-invariant R
        isRTimeInvariant = true;
        % Output estimated model state x
        showEstimatedState = true;
        % Output estimated input
        showEstimatedInputs = true;
        % Output estimated model output y
        showEstimatedOutput = false;
        % Output state estimation covariance P
        showStateCovariance = false;
    end

    % Public, nontunable properties
    properties(Nontunable)
        % A
        AMat = 1;
        % B
        BMat = 1;
        % H
        HMat = 1;
        % C
        CMat = 1;
        % D
        DMat = 1;
        % G
        GMat = 1;
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
        % Model source
        ModelSourceChoice = 'Individual matrices';
    end
    
    % Public, tunable, positive integre properties
    properties(PositiveInteger)
        % Number of states
        NxUserDefined = 1;
        % Number of measurements
        NyUserDefined = 1;
        % Number of known inputs
        NuUserDefined = 1;
        % Number of unknown inputs
        NeUserDefined = 1;
    end
    
    % Discrete states
    properties(DiscreteState)
        xHat;
        PHat;
        prevInput;
    end
    
    % Hidden constants
    properties(Hidden, Constant)
        ModelSourceChoiceSet = matlab.system.StringSet({...
            'Individual matrices',...
            'Input port'});
    end

    % Pre-computed positive-integer constants
    properties(Access = private, PositiveInteger)
        Nx;
        Ny;
        Nu;
        Ne;
    end

    methods(Access = protected)
        %% setupImpl & stepImpl & resetImpl
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            
            % Define size of the system
            [obj.Nx, obj.Nu, obj.Ne, obj.Ny] = getSystemSize(obj);
        end

        function varargout = stepImpl(obj, varargin)
            % Rename input
            num = 0;
            % Add input for known inputs
            if obj.hasKnownInput
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
            % Add input for state matrices
            if strcmp(obj.ModelSourceChoice, 'Input port')
                if obj.hasKnownInput
                    num = num + 6;  % Matrices A, B, H, C, D, and G
                    A = varargin{num - 5};
                    B = varargin{num - 4};
                    H = varargin{num - 3};
                    C = varargin{num - 2};
                    D = varargin{num - 1};
                    G = varargin{num};
                else
                    num = num + 4;  % Matrices A, H, C, and G
                    A = varargin{num - 3};
                    B = zeros(size(obj.Nx, 1), 1);
                    H = varargin{num - 2};
                    C = varargin{num - 1};
                    D = zeros(size(obj.Ny, 1), 1);
                    G = varargin{num};
                end
            else
                if obj.hasKnownInput
                    A = obj.AMat;
                    B = obj.BMat;
                    H = obj.HMat;
                    C = obj.CMat;
                    D = obj.DMat;
                    G = obj.GMat;
                else
                    A = obj.AMat;
                    B = zeros(size(obj.Nx, 1), 1);
                    H = obj.HMat;
                    C = obj.CMat;
                    D = zeros(size(obj.Ny, 1), 1);
                    G = obj.GMat;
                end
            end
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
            [xCurrGPrev, PCurrGPrev] = obj.aPrioriUpdate(A, B, Q, ...
                xPrevGPrev, PPrevGPrev, uPrev);
            
            % A posteriori update
            [xCurrGCurr, PCurrGCurr, L] = obj.aPosterioriUpdate(C, D, ...
                G, H, R, xCurrGPrev, PCurrGPrev, uCurr, yCurr);
            
            % Compute unknown input (if necessary)
            if obj.showEstimatedInputs || obj.showEstimatedOutput
                eCurr = obj.estimateInput(C, D, H, G, L, ...
                    xCurrGPrev, xCurrGCurr, uCurr, yCurr);
            end
            
            % Compute estimate of the output (if necessary)
            if obj.showEstimatedOutput
                yCurrGCurr = obj.estimateOutput(C, D, G, ...
                    xCurrGCurr, uCurr, eCurr);
            end
            
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
            % Output estimated input
            if obj.showEstimatedInputs
                varargout = [varargout{:}, {eCurr}];
            end
            % Output estimated model output y
            if obj.showEstimatedOutput
                varargout = [varargout{:}, {yCurrGCurr}];
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
            
            % Check matrices size if matrices are provided as properties
            if strcmp(obj.ModelSourceChoice, 'Individual matrices')
                obj.checkSystemMatricesSize();
            end
            
            % Check sampling time property
            if ~(obj.Ts == -1 || obj.Ts > 0)
                error('The sampling time must be positive or -1')
            end
        end
    end
    
    methods(Access = private)
        function checkSystemMatricesSize(obj)
            % Throw an error if the size of the system matrices are
            % incompatible
            
            errorMsg = '';
            
            % A is square
            if size(obj.AMat,1) ~= size(obj.AMat,2)
                errorMsg = [errorMsg 'The A matrix must be square. '];
            end
            
            % Check number of rows
            if obj.hasKnownInput && size(obj.BMat,1) ~= size(obj.AMat,1)
                errorMsg = [errorMsg ...
                    'The A and B matrices must have the same number '...
                    'of rows. '];
            end
            if size(obj.HMat,1) ~= size(obj.AMat,1)
                errorMsg = [errorMsg ...
                    'The A and H matrices must have the same number '...
                    'of rows. '];
            end
            if obj.hasKnownInput && size(obj.DMat,1) ~= size(obj.CMat,1)
                errorMsg = [errorMsg ...
                    'The C and D matrices must have the same number '...
                    'of rows. '];
            end
            if size(obj.GMat,1) ~= size(obj.CMat,1)
                errorMsg = [errorMsg ...
                    'The C and G matrices must have the same number '...
                    'of rows. '];
            end
            
            % Check number of columns
            if size(obj.AMat,2) ~= size(obj.CMat,2)
                errorMsg = [errorMsg ...
                    'The A and C matrices must have the same number '...
                    'of columns. '];
            end
            if obj.hasKnownInput && size(obj.BMat,2) ~= size(obj.DMat,2)
                errorMsg = [errorMsg ...
                    'The B and D matrices must have the same number '...
                    'of columns. '];
            end
            if size(obj.HMat,2) ~= size(obj.GMat,2)
                errorMsg = [errorMsg ...
                    'The H and G matrices must have the same number '...
                    'of columns. '];
            end
            
            % Throw error if necessary
            if ~isempty(errorMsg)
                error(errorMsg)
            end
        end
        
        function [Nx, Nu, Ne, Ny] = getSystemSize(obj)
            if strcmp(obj.ModelSourceChoice, 'Input port')
                Nx = obj.NxUserDefined;
                Ny = obj.NyUserDefined;
                if obj.hasKnownInput
                    Nu = obj.NuUserDefined;
                else
                    Nu = 1;
                end
                Ne = obj.NeUserDefined;
            elseif strcmp(obj.ModelSourceChoice, 'Individual matrices')
                Nx = size(obj.AMat, 1);
                Ny = size(obj.CMat, 1);
                if obj.hasKnownInput
                    Nu = size(obj.BMat, 2);
                else
                    Nu = 1;
                end
                Ne = size(obj.HMat, 2);
            end
        end
        
        function [xCurrGPrev, PCurrGPrev] = aPrioriUpdate(~, A, B, Q,...
                xPrevGPrev, PPrevGPrev, uPrev)
            % A priori update
            xCurrGPrev = A * xPrevGPrev + B * uPrev;
            PCurrGPrev = A * PPrevGPrev * A' + Q;
        end
        
        function [xCurrGCurr, PCurrGCurr, L] = aPosterioriUpdate(obj, C,...
                D, G, H, R, xCurrGPrev, PCurrGPrev, uCurr, yCurr)
            % Filter update
            Rtilde = C * PCurrGPrev * C' + R;
            F = PCurrGPrev * C';
            if prod(H == 0,'all') && prod(G == 0,'all') % H  = 0 and G  = 0
                L = F / Rtilde;
            elseif prod(H == 0,'all')                   % H  = 0 and G ~= 0
                Phi = G;
                Ome = -F / Rtilde * G;
                L = (F + Ome / (Phi' / Rtilde * Phi) * Phi') / Rtilde;
            elseif prod(G == 0,'all')                   % H ~= 0 and G  = 0
                V = C*H;
                Pi = (V' / Rtilde * V) \ V' / Rtilde;
                L = H*Pi + F*Rtilde*(eye(obj.Ny) - V*Pi);
            else                                        % H ~= 0 and G ~= 0
                % Find non-zero columns of H and G
                nzH_col = find(sum(abs(H)) ~= 0); % Non-zero columns of H
                nzG_col = find(sum(abs(G)) ~= 0); % Non-zero columns of G
                NnzG = numel(nzG_col);    % Number of non-zero columns of G
                % Compute filter
                Phi = [-G(:,nzG_col) C*H(:,nzH_col)];
                Ome = [zeros(obj.Nx, NnzG) H(:,nzH_col)] - F/Rtilde*Phi;
%                 L = (F + Ome / (Phi' / Rtilde * Phi) * Phi') / Rtilde;
                L = (F + Ome * pinv(Phi' / Rtilde * Phi) * Phi') / Rtilde; % This seem to avoid badly scaed matrix warning
            end
            
            % A posteriori update
            xCurrGCurr = xCurrGPrev + L * (...
                yCurr - (C * xCurrGPrev + D * uCurr) ...
            );
            PCurrGCurr = L * Rtilde * L' - F * L' - L * F' + PCurrGPrev;
        end
        
        function eCurr = estimateInput(~, C, D, H, G, L, ...
                xCurrGPrev, xCurrGCurr, uCurr, yCurr)
            % Diagonalize pinv(H)*H and pinv(G)*G
            [VH, SigmaH] = eig(pinv(H)*H);
            [VG, SigmaG] = eig(pinv(G)*G);
            VH = inv(VH);
            VG = inv(VG);
            % Create vector of eigenvalues of pinv(H)*H and pinv(G)*G
            eigH = diag(SigmaH);
            eigG = diag(SigmaG);
            % Find non-zero eigenvalues
            nzEigH = find(abs(eigH) >= eps);
            nzEigG = find(abs(eigG) >= eps);
            % Compute the unknown inputs using the two equations
            eH = pinv(H) * L * (yCurr - (C * xCurrGPrev + D * uCurr));
            eG = pinv(G) *     (yCurr - (C * xCurrGCurr + D * uCurr));
            % Solve the system to compute the unknown inputs
            eCurr = pinv([VH(nzEigH,:); VG(nzEigG,:)]) * ...
                [VH(nzEigH,:) * eH; VG(nzEigG,:) * eG];
        end
        
        function yCurrGCurr = estimateOutput(~, C, D, G, ...
                xCurrGCurr, uCurr, eCurrGCurr)
            % Compute output
            yCurrGCurr = C*xCurrGCurr + D*uCurr + G*eCurrGCurr;
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
            if obj.hasKnownInput
                num = num + 1;
            end
            % Add input for measurement
            num = num + 1;
            % Add input for state matrices
            if strcmp(obj.ModelSourceChoice, 'Input port')
                if obj.hasKnownInput
                    num = num + 6;  % Matrices A, B, H, C, D, and G
                else
                    num = num + 4;  % Matrices A, H, C, and G
                end
            end
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
            % Output estimated input
            if obj.showEstimatedInputs
                num = num + 1;
            end
            % Output estimated model output y
            if obj.showEstimatedOutput
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
            if obj.hasKnownInput
                varargout = [varargout{:}, {'u'}];
            end
            % Add input for measurement
            varargout = [varargout{:}, {'y'}];
            % Add input for state matrices
            if strcmp(obj.ModelSourceChoice, 'Input port')
                if obj.hasKnownInput
                    % Matrices A, B, H, C, D, and G
                    varargout = [varargout{:}, ...
                        {'A'}, {'B'}, {'H'}, {'C'}, {'D'}, {'G'}];
                else
                    % Matrices A, H, C, and G
                    varargout = [varargout{:}, {'A'}, {'H'}, {'C'}, {'G'}];
                end
            end
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
            % Output estimated input
            if obj.showEstimatedInputs
                varargout = [varargout{:}, {'ehat'}];
            end
            % Output estimated model output y
            if obj.showEstimatedOutput
                varargout = [varargout{:}, {'yhat'}];
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
            % Output estimated input e
            if obj.showEstimatedInputs
                varargout = [varargout{:}, {true}];
            end
            % Output estimated model output y
            if obj.showEstimatedOutput
                varargout = [varargout{:}, {true}];
            end
            % Output state estimation covariance P
            if obj.showStateCovariance
                varargout = [varargout{:}, {true}];
            end
        end

        function varargout = getOutputSizeImpl(obj)
            % Return size for each output port
            
            [Nbx, ~, Nbe, Nby] = getSystemSize(obj);
            
            varargout = {};
            % Output estimated model state x
            if obj.showEstimatedState
                varargout = [varargout{:}, {[Nbx, 1]}];
            end
            % Output estimated input
            if obj.showEstimatedInputs
                varargout = [varargout{:}, {[Nbe, 1]}];
            end
            % Output estimated model output y
            if obj.showEstimatedOutput
                varargout = [varargout{:}, {[Nby, 1]}];
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
            % Output estimated input
            if obj.showEstimatedInputs
                varargout = [varargout{:}, {'double'}];
            end
            % Output estimated model output y
            if obj.showEstimatedOutput
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
            % Output estimated input
            if obj.showEstimatedInputs
                varargout = [varargout{:}, {false}];
            end
            % Output estimated model output y
            if obj.showEstimatedOutput
                varargout = [varargout{:}, {false}];
            end
            % Output state estimation covariance P
            if obj.showStateCovariance
                varargout = [varargout{:}, {false}];
            end
        end
        
        %% Set inactive properties
        function flag = isInactivePropertyImpl(obj,propertyName)
            
            % State matrices are required only if the option 
            % indivudual matrices are used
            if ismember(propertyName, {'AMat','CMat','HMat','GMat'})
                flag = not(strcmp(obj.ModelSourceChoice,...
                    'Individual matrices'));
            % Only if individual matrices are provided and there is an
            % input
            elseif ismember(propertyName, {'BMat','DMat'})
                flag = not(strcmp(obj.ModelSourceChoice,...
                    'Individual matrices') && obj.hasKnownInput);
            % Sizes are required only if the option input port is used
            elseif ismember(propertyName, {'NxUserDefined',...
                    'NeUserDefined','NyUserDefined'})
                flag = not(strcmp(obj.ModelSourceChoice,...
                    'Input port'));
            % Only if input port and there is an input
            elseif ismember(propertyName, {'NuUserDefined'})
                flag = not(strcmp(obj.ModelSourceChoice,...
                    'Input port') && obj.hasKnownInput);
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
            icon = {'Unbiased','Minimum','Variance','Filter'};
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
            header = matlab.system.display.Header('UMV',...
                'Title','Unbiased Minimum Variance Filter',...
                'Text', 'Estimates the states and unknown inputs of a discrete-time system.');
        end

        function groups = getPropertyGroupsImpl
            % Customize dialog box
            
            % +-----------------------+
            % |    Model Parameters   |
            % +-----------------------+
            % Create section for system model
            SystemModelSection = matlab.system.display.Section(...
                'Title','System model',...
                'PropertyList',{'ModelSourceChoice',...
                'AMat','BMat','HMat','CMat','DMat','GMat',...
                'NxUserDefined','NuUserDefined','NeUserDefined',...
                'NyUserDefined',...
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
                'PropertyList',{'hasKnownInput'});
            
            % Create section for optional output ports
            additionalOutportsSection = matlab.system.display.Section(...
                'Title','Additional Outports', ...
                'PropertyList',{'showEstimatedState',...
                'showEstimatedInputs','showEstimatedOutput',...
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
