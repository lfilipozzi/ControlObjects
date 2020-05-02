classdef LtvMpc < matlab.System & ...
    matlab.system.mixin.CustomIcon & ...
    matlab.system.mixin.Propagates & ...
    matlab.system.mixin.SampleTime
% LTVMPC provides an implementation of an augmented LTV-MPC controller for 
% reference tracking in Simulink.
% 
% The augmented LTV-MPC problem is defined as follow:
%   Minimize:
%       - Reference tracking error over prediction horizon.
%       - Perturbed control input over control horizon.
%       - Absolute control input over control horizon (this term is not
%       present in the standard formulation of the augmented LTV-MPC).
%       - Slack variable penalty.
%   Subject to:
%       - Linearized state-space difference equation (linearized at the
%       previous absolute control input and current state). The plant can
%       include measured disturbances (they should correspond to the last
%       control inputs of the state-space model).
%       - Constraints on the perturbed state and perturbed input (possibly
%       soft constraints).
%       - Bounds on the absolute control input.
% 
% This MATLAB System Object support the following features:
%  - Measured disturbance in the MPC model
%  - Soft constraints
%  - Different control and prediction horizon
%  - Use of estimated signals in the MPC problem formulation
%  - Scaling factors for the optimization manipulated variables (control 
%    inputs and slack variables) to scale the QP problem and enhance the 
%    solver performance. The same scaling factor is used for all the slack 
%    variables.
% 
% The MPC formulation requires mainly three functions to define the MPC
% problem.
%  - A function returning the jacobian matrices in the continuous-time
%    domain. The prototype of the function if as follow
%       [Ac, Bc, Cc, Dc, yop] = getJacobianMat(xop, uop, estim)
%    where Ac, Bc, Cc, and Dc are the continuous state-space matrices of
%    the MPC model; xop, uop, and yop are the state, input and output 
%    operating points; and estim is a structure containing estimated 
%    signals.
%  - A function returning the stage cost. Its prototype is
%       [Qe, Rd, Ru] = obj.getCostFuncMat(xop, uop, estim);
%    where Qe, Rd, and Ru defines the stage cost; xop and uop are 
%    the state and input operating points; and estim is a structure 
%    containing estimated signals.
%    Remark: The stage cost is defined as follow
%         J_k = e_k^T Qe e_k + Delta u_k^T Rd Delta u_k + u_k^t Ru u_k
%    where e_k is the reference tracking error, u_k is the control input
%    and Delta u_k is the perturbed control input.
%  - If the MPC has inequlity constraints, the matrices defining the
%    polytopic constraints are given by the following function
%       [Ax, Au, bop] = getPolytopeMat(xop, uop, estim)
%    where Ax, and Au deines the inequality constraints; bop is the value 
%    of the constraints evaluated at (xop, uop); xop and uop are the state
%    and input operating points and estim is a structure containing 
%    estimated signals.
% 
% @author Louis Filipozzi

properties(Nontunable,Logical)
    % Show optimal cost
    showCost = false;
    % Show solving time
    showSolveTime = false;
    % Show optimal control sequence
    showControlSequence = false;
    % Show predicted state sequence
    showStateSequence = false;
    % Show predicted output sequence
    showOutputSequence = false;
    % Show slack variable
    showSlackVariable = false;
    % Measured disturbance
    hasMeasuredDisturbance = false;
    % Inequality constrainsts
    hasIneqConstraints = false;
    % Actuator bounds
    hasActuatorBounds = false;
    % Estimated signals
    hasEstimatedSignals = false;
    % Soft inequality constraints
    hasSoftConstraints = false;
    % Add scaling factors
    hasScaling = false;
    % Discrete model
    isModelDiscrete = false;
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
    % Number of reference signals
    Nr = 2;
end

properties(Nontunable)
    % Number of measured disturbance
    Nw = 0;
    % Sampling time
    Ts = 0.01;
    % Cost function
    getCostFuncMat = @(xop, uop, estim) 0;
    % Jacobian function
    getJacobianMat = @(xop, uop, estim) 1;
    % Inequality matrices function
    getPolytopeMat = @(xop, uop, estim) {[],[]};
    % Indices
    softConstraintsIndex = {0};
    % Weight
    softConstraintsWeight = 1;
    % Control inputs and measured disturbances
    actuatorScaling = 1;
    % Slack variables
    slackScaling = 1;
    % QP solver
    qpSolverChoice = 'quadprog';
end

properties(Hidden, Constant)
    qpSolverChoiceSet = matlab.system.StringSet({...
        'quadprog','qpOASES','qpOASES_sequence'...
    });
end

properties(Access = private)
    mpcController LtvMpcController;
    qpSolver IQpSolver = QuadprogSolver();
end

methods(Access = protected)
    function setupImpl(obj)
        % Perform one-time calculations, such as computing constants
        
        % If no measured disturbance, set Nw to 0 automatically
        if ~obj.hasMeasuredDisturbance
            obj.Nw = 0;
        end

        % Set QP solver options
        if strcmp(obj.qpSolverChoice, 'quadprog')
            obj.qpSolver = QuadprogSolver();
        elseif strcmp(obj.qpSolverChoice, 'qpOASES')
            obj.qpSolver = QpOasesSolver();
        elseif strcmp(obj.qpSolverChoice, 'qpOASES_sequence')
            obj.qpSolver = QpOasesSequenceSolver();
        end
        
        % Create MPC controller
        obj.mpcController = LtvMpcController(...
            obj.Nt,obj.Np,obj.Nx,obj.Nu,obj.Nr,obj.Nw,obj.Ts,...
            obj.qpSolver...
        );
    
        if obj.isModelDiscrete
            obj.mpcController.setPlantModelSamplingTime(obj.Ts);
        else
            obj.mpcController.setPlantModelSamplingTime(0);
        end
        
        if obj.hasIneqConstraints && obj.hasSoftConstraints
            obj.mpcController.setSoftConstraints(...
                obj.softConstraintsIndex, obj.softConstraintsWeight...
        );
        else
            obj.hasSoftConstraints = false;
            obj.showSlackVariable = false;
            obj.mpcController.setSoftConstraints({},[]);
        end
        
        if obj.hasScaling
            obj.mpcController.setScaleFactors(...
                obj.actuatorScaling, obj.slackScaling...
        );
        end
        
        obj.mpcController.setup();
    end

    function varargout = stepImpl(obj,varargin)
        % Implement algorithm.

        % Rename inputs
        ref    = varargin{1};
        states = varargin{2};
        num = 2;
        if obj.hasEstimatedSignals
            num = num + 1;
            estim = varargin{num};
        else
            estim = {};
        end
        if obj.hasMeasuredDisturbance
            num = num + 1;
            wMeas = varargin{num};
        else
            wMeas = [];
        end
        if obj.hasActuatorBounds
            num = num + 2;
            umin = varargin{num - 1};
            umax = varargin{num};
        else
            umin = [];
            umax = [];
        end
        if obj.hasIneqConstraints
            num = num + 1;
            b = varargin{num};
        else
            b = [];
        end

        % Setup chronometer
        if obj.showSolveTime; tic; end
        xop = states;
        uPrev = obj.mpcController.getPreviousInput();
        uop = [uPrev; wMeas];
        [Qe,Rd,Ru] = obj.getCostFuncMat(xop, uop, estim);
        [A,B,C,D,yop]  = obj.getJacobianMat(xop, uop, estim);
        if (obj.hasIneqConstraints)
            [Ax,Au,bop]    = obj.getPolytopeMat(xop, uop, estim);
        else
            Ax  = [];
            Au  = [];
            bop = [];
        end
        obj.mpcController.setCostFunction(Qe,Rd,Ru);
        obj.mpcController.setMeasuredDisturbance(wMeas);
        obj.mpcController.setReferenceTarget(ref);
        obj.mpcController.setPlantModel(A,B,C,D,xop,yop);
        obj.mpcController.setConstraints(Ax,Au,bop,b);
        obj.mpcController.setActuatorBounds(umin,umax);

        % Solve MPC problem using batch approach
        [u_seq, x_seq, y_seq, cost, exitflag, slack] = obj.mpcController.step();

        % Warning if the problem has not been solved
        errMessage = obj.qpSolver.getWarningMessage(exitflag);
        if ~isempty(errMessage)
            errMessage = ['Solver failed at time ',...
                num2str(obj.getCurrentTime(),10),'! ',errMessage];
            warning(errMessage);
        end

        % Extract first control input
        u = u_seq(1:obj.Nu,1);

        % Record time to solve
        if obj.showSolveTime
            solveTime = toc;
        end

        % Returns output: [u, solveTime, cost, exitflag]
        num = 1;
        varargout{1} = u;
        if obj.showSolveTime
            num = num + 1;
            varargout{num} = solveTime;
        end
        if obj.showCost
            num = num + 1;
            varargout{num} = cost;
        end
        if obj.showControlSequence
            num = num + 1;
            varargout{num} = u_seq;
        end
        if obj.showStateSequence
            num = num + 1;
            varargout{num} = x_seq;
        end
        if obj.showOutputSequence
            num = num + 1;
            varargout{num} = y_seq;
        end
        num = num + 1;
        varargout{num} = exitflag;
        if obj.hasIneqConstraints && obj.hasSoftConstraints && ...
                obj.showSlackVariable
            num = num + 1;
            varargout{num} = slack;
        end
    end

    function resetImpl(obj)
        % Initialize / reset discrete-state properties
        obj.mpcController.reset();
    end

    function releaseImpl(obj)
        % Release resources
        obj.mpcController.release();
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

methods(Access = protected)
    function sts = getSampleTimeImpl(obj)
        % Define sampling time
        sts = createSampleTime(obj,'Type','Discrete',...
          'SampleTime',obj.Ts,'OffsetTime',0);
    end

    % Indicate number of input ports of the system object
    function num = getNumInputsImpl(obj)
        % Two inputs minimum for state feedback and reference signal
        num = 2;
        if obj.hasEstimatedSignals
            % Add inport for estimated signals
            num = num + 1;
        end
        if obj.hasMeasuredDisturbance
            % Add inport for measured disturbance
            num = num + 1;
        end
        if obj.hasActuatorBounds
            % Add two signals for max and min actuator bounds
            num = num + 2;
        end
        if obj.hasIneqConstraints
            % Add one input for ineqaulity bounds
            num = num + 1;
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
        if obj.showStateSequence
            num = num + 1;
        end
        if obj.showOutputSequence
            num = num + 1;
        end
        if obj.hasIneqConstraints && obj.hasSoftConstraints && ...
                obj.showSlackVariable
            num = num + 1;
        end
    end
    
    function varargout = getInputNamesImpl(obj)
        % Non-optional inputs
        num = 1;
        varargout{num} = 'reference';
        num = num + 1;
        varargout{num} = 'state';
        % Optional inputs
        if obj.hasEstimatedSignals
            num = num + 1;
            varargout{num} = 'estim. signals';
        end
        if obj.hasMeasuredDisturbance
            num = num + 1;
            varargout{num} = 'meas. dist.';
        end
        if obj.hasActuatorBounds
            num = num + 2;
            varargout{num-1} = 'min. actuator';
            varargout{num}   = 'max. actuator';
        end
        if obj.hasIneqConstraints
            num = num + 1;
            varargout{num}= 'ineq. bounds';
        end
    end
    
    function varargout = getOutputNamesImpl(obj)
        % Non-optional outputs
        num = 1;
        varargout{num} = 'control';
        % Optional outputs
        if obj.showSolveTime
            num = num + 1;
            varargout{num} = 'time';
        end
        if obj.showCost
            num = num + 1;
            varargout{num} = 'cost';
        end
        if obj.showControlSequence
            num = num + 1;
            varargout{num} = 'controls seq.';
        end
        if obj.showStateSequence
            num = num + 1;
            varargout{num} = 'pred. states';
        end
        if obj.showOutputSequence
            num = num + 1;
            varargout{num} = 'pred. outputs';
        end
        num = num + 1;
        varargout{num} = 'exitflag';
        if obj.hasIneqConstraints && obj.hasSoftConstraints && ...
                obj.showSlackVariable
            num = num + 1;
            varargout{num} = 'slack var.';
        end
    end
    
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
        num = 1;
        varargout{num} = true;

        % solving time
        if obj.showSolveTime
            num = num + 1;
            varargout{num} = true;
        end

        % optimal cost
        if obj.showCost
            num = num + 1;
            varargout{num} = true;
        end

        % control sequence
        if obj.showControlSequence
            num = num + 1;
            varargout{num} = true;
        end
        
        % predicted state sequence
        if obj.showStateSequence
            num = num + 1;
            varargout{num} = true;
        end
        
        % predicted output sequence
        if obj.showOutputSequence
            num = num + 1;
            varargout{num} = true;
        end

        % exitflag
        num = num + 1;
        varargout{num} = true;

        % slack variable
        if obj.hasIneqConstraints && obj.hasSoftConstraints && ...
                obj.showSlackVariable
            num = num + 1;
            varargout{num} = true;
        end
    end

    function varargout = getOutputSizeImpl(obj)
        % Return size for each output port

        % control
        num = 1;
        varargout{num} = [obj.Nu 1];

        % solving time
        if obj.showSolveTime
            num = num + 1;
            varargout{num} = [1 1];
        end

        % optimal cost
        if obj.showCost
            num = num + 1;
            varargout{num} = [1 1];
        end

        % control sequence
        if obj.showControlSequence
            num = num + 1;
            varargout{num} = [obj.Nu obj.Nt];
        end
        
        % predicted state sequence
        if obj.showStateSequence
            num = num + 1;
            varargout{num} = [obj.Nx obj.Np];
        end
        
        % predicted output sequence
        if obj.showOutputSequence
            num = num + 1;
            varargout{num} = [obj.Nr obj.Np];
        end

        % exitflag
        num = num + 1;
        varargout{num} = [1 1];

        % slack variables
        if obj.hasIneqConstraints && obj.hasSoftConstraints && ...
                obj.showSlackVariable
            num = num + 1;
            varargout{num} = [numel(obj.softConstraintsIndex) obj.Np];
        end
    end

    function varargout = getOutputDataTypeImpl(obj)
        % u_output
        num = 1;
        varargout{num} = 'double';

        % solving time
        if obj.showSolveTime
            num = num + 1;
            varargout{num} = 'double';
        end

        % optimal cost
        if obj.showCost
            num = num + 1;
            varargout{num} = 'double';
        end

        % control sequence
        if obj.showControlSequence
            num = num + 1;
            varargout{num} = 'double';
        end
        
        % predicted state sequence
        if obj.showStateSequence
            num = num + 1;
            varargout{num} = 'double';
        end
        
        % predicted output sequence
        if obj.showOutputSequence
            num = num + 1;
            varargout{num} = 'double';
        end

        % exitflag
        num = num + 1;
        varargout{num} = 'double';

        % slack variables
        if obj.hasIneqConstraints && obj.hasSoftConstraints && ...
                obj.showSlackVariable
            num = num + 1;
            varargout{num} = 'double';
        end
    end

    function varargout = isOutputComplexImpl(obj)
        % u_output
        num = 1;
        varargout{num} = false;

        % solving time
        if obj.showSolveTime
            num =  num + 1;
            varargout{num} = false;
        end

        % optimal cost
        if obj.showCost
            num = num + 1;
            varargout{num} = false;
        end

        % control sequence
        if obj.showControlSequence
            num = num + 1;
            varargout{num} = false;
        end
        
        % predicted state sequence
        if obj.showStateSequence
            num = num + 1;
            varargout{num} = false;
        end
        
        % predicted output sequence
        if obj.showOutputSequence
            num = num + 1;
            varargout{num} = false;
        end

        % exitflag
        num = num + 1;
        varargout{num} = false;

        % slack variable
        if obj.hasIneqConstraints && obj.hasSoftConstraints && ...
                obj.showSlackVariable
            num = num + 1;
            varargout{num} = false;
        end
    end
    
    function flag = isInactivePropertyImpl(obj,propertyName)
        % getPolytopeMat is useful only if we have inequality 
        % constraints
        if ismember(propertyName, {'getPolytopeMat'})
            flag = not(obj.hasIneqConstraints);
        % Nw is needed only if we have measured disturbance
        elseif ismember(propertyName, {'Nw'})
            flag = not(obj.hasMeasuredDisturbance);
        elseif ismember(propertyName, {'hasSoftConstraints'})
            flag = not(obj.hasIneqConstraints);
        % Soft constraints properties only with soft constraints
        elseif ismember(propertyName, {'softConstraintsIndex',...
                'softConstraintsWeight','showSlackVariable'})
            flag = not(obj.hasIneqConstraints && obj.hasSoftConstraints);
        % Scaling used for QP solver
        elseif ismember(propertyName, {'actuatorScaling'})
            flag = not(obj.hasScaling);
        elseif ismember(propertyName, {'slackScaling'})
            flag = not(obj.hasScaling && obj.hasSoftConstraints && ...
                obj.hasIneqConstraints);
        % All other properties are active
        else
            flag = false;
        end
    end
    
    function icon = getIconImpl(~)
        icon = {'Linear Time','Varying','MPC'};
    end
end 
    
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
        header = matlab.system.display.Header('LtvMpc',...
            'Title','Linear Time-Varying MPC',...
            'Text', ['MATLAB System implementation of an LTV-MPC. ',...
            'The MPC is designed for reference tracking and uses ',...
            'an augmented formulation. Use MATLAB function to ',...
            'provide the cost function matrices, jacobian matrices ',...
            'and inequality matrices.']);
    end

    function groups = getPropertyGroupsImpl
        % Customize dialog box

        % +-----------------------+
        % |    MPC Formulation    |
        % +-----------------------+
        % Create section for model
        MPCModelSection = matlab.system.display.Section(...
            'Title','Model',...
            'PropertyList',{'getJacobianMat','isModelDiscrete'});

        % Create section for cost function
        MPCCostSection = matlab.system.display.Section(...
            'Title','Cost',...
            'PropertyList',{'getCostFuncMat'});

        MPCConstraints = matlab.system.display.Section(...
            'Title','Constraints',...
            'PropertyList',{'hasActuatorBounds',...
            'hasIneqConstraints','getPolytopeMat',...
            'hasSoftConstraints','softConstraintsIndex',...
            'softConstraintsWeight'});

        % Create tab for MPC formulation
        MPCFormulationGroup = matlab.system.display.SectionGroup(...
            'Title','Formulation',...
            'Sections',[MPCModelSection, MPCCostSection, ...
            MPCConstraints]);

        % +-----------------------+
        % |    MPC parameters     |
        % +-----------------------+
        % Crete sectio for QP solver
        qpSolverSection = matlab.system.display.Section(...
            'Title','QP Solver',...
            'PropertyList',{'qpSolverChoice'});

        % Create section for hyper-parameters
        hyperParametersSection = matlab.system.display.Section(...
            'Title','Hyper-parameters', ...
            'PropertyList',{'Ts','Nt','Np'});

        % Create section for sizes
        sizesSection = matlab.system.display.Section(...
            'Title','Signals Sizes', ...
            'PropertyList',{'Nx','Nu','Nr'});

        % Create section for scaling the QP problem
        scalingSection = matlab.system.display.Section(...
            'Title','Scaling the Quadratic Problem', ...
            'PropertyList',{'hasScaling','actuatorScaling',...
            'slackScaling'});

        % Create tab for MPC parameters configuration
        MPCParametersGroup = matlab.system.display.SectionGroup(...
            'Title','Parameters',...
            'Sections',[qpSolverSection, hyperParametersSection, ...
            sizesSection, scalingSection]);

        % +-----------------------+
        % |   MPC optional ports  |
        % +-----------------------+
        % Create section for optional input ports
        additionalInportsSection = matlab.system.display.Section(...
            'Title','Additional Inports', ...
            'PropertyList',{'hasMeasuredDisturbance','Nw',...
            'hasEstimatedSignals'});

        % Create section for optional output ports
        additionalOutportsSection = matlab.system.display.Section(...
            'Title','Additional Outports', ...
            'PropertyList',{'showCost','showSolveTime',...
            'showControlSequence','showStateSequence',...
            'showOutputSequence','showSlackVariable'});

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
