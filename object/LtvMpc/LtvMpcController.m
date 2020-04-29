classdef LtvMpcController < handle
% LTVMPCCONTROLLER Define an augmented LTV-MPC controller for reference
% tracking.
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

properties(Access = private)
    % Control horizon
    Nt = 10;
    % Prediction horizon
    Np = 10;
    % Number of states
    Nx = 1;
    % Number of inputs
    Nu = 1;
    % Number of reference signals
    Nr = 1;
    % Number of measured disturbance
    Nw = 0;
    % Sampling time
    Ts = 0.01;
    % Cost function matrices
    cost = struct('Qe',1,'Rd',1,'Ru',0);
    % Continuous state-space matrices
    plant = struct('A',1,'B',1,'C',1,'D',1,'yop',1,'Ts',0);
    % Control input at the last timestep
    uPrev = 0;
    % Measured disturbance
    wMeas = [];
    % Reference target
    yRef = 0;
    % Inequality contraints
    constraints = struct('Ax',[],'Au',[],'bop',[],'b',[]);
    % Soft constraints
    softConstraints = struct('index',{},'weight',[]);
    % Scale factors
    scaling = struct('input',[],'slack',[],'factors',[]);
    % Input bounds
    bounds = struct('lb',[],'ub',[]);
    % QP solver
    qpSolver IQpSolver = QuadprogSolver();
end

methods(Access = public)
    function obj = LtvMpcController(Nt, Np, Nx, Nu, Nr, Nw, Ts, qpSolver)
        % Contructor of the class
        
        % Run checks
        if ~isa(qpSolver,'IQpSolver')
            error('qpSolver must implement the IQpSolver interface');
        end
        
        % Set attributes
        obj.Nt = Nt;
        obj.Np = Np;
        obj.Nx = Nx;
        obj.Nu = Nu;
        obj.Nr = Nr;
        obj.Nw = Nw;
        obj.Ts = Ts;
        obj.qpSolver = qpSolver;
    end
    
    function setCostFunction(obj,Qe,Rd,Ru)
        % Set the cost function of the MPC problem
        obj.cost.Qe = Qe;
        obj.cost.Rd = Rd;
        obj.cost.Ru = Ru;
    end
    
    function setMeasuredDisturbance(obj, w)
        % Set the measured disturbance of the MPC problem
        obj.wMeas = w;
    end
    
    function setReferenceTarget(obj,y)
        % Set the target reference of the MPC problem
        obj.yRef = y;
    end
    
    function setPlantModel(obj,A,B,C,D,yop)
        % Set the continuous state-space model
        obj.plant.A   = A;
        obj.plant.B   = B;
        obj.plant.C   = C;
        obj.plant.D   = D;
        obj.plant.yop = yop;
    end
    
    function setPlantModelSamplingTime(obj,Ts)
        obj.plant.Ts = Ts;
    end
    
    function setConstraints(obj,Ax,Au,bop,b)
        % Set the polytopic constraints of the MPC problem
        if ~isempty(Ax) && ~isempty(Au) && ~isempty(b)
            obj.constraints.Ax  = Ax;
            obj.constraints.Au  = Au;
            obj.constraints.bop = bop;
            obj.constraints.b   = b;
        else
            obj.constraints.Ax  = zeros(0,obj.Nx);
            obj.constraints.Au  = zeros(0,obj.Nu+obj.Nw);
            obj.constraints.bop = [];
            obj.constraints.b   = [];
        end
    end
    
    function setActuatorBounds(obj,lb,ub)
        % Set the actuator bounds of the MPC problem
        if ~isempty(lb)
            obj.bounds.lb = lb;
        else
            obj.bounds.lb = -inf(obj.Nu,1);
        end
        if ~isempty(ub)
            obj.bounds.ub = ub;
        else
            obj.bounds.ub = inf(obj.Nu,1);
        end
    end
    
    function setSoftConstraints(obj,index,weight)
        % Specify soft constraints
        obj.softConstraints = struct('index',{index},'weight',weight);
    end
    
    function setScaleFactors(obj,input,slack)
        % Set the scaling factors
        obj.scaling = struct('input',input,'slack',slack);
    end
    
    function uPrev = getPreviousInput(obj)
        % Return the previous control input
        uPrev = obj.uPrev;
    end
    
    function setup(obj)
        % Setup the MPC controller (perform one-time calculations)
        obj.qpSolver.setup();
    end
    
    function [u_seq_opt, cost, exitflag, slack] = step(obj)
        % Run one optmization of the MPC controller. Return the optimal
        % control sequence, its cost, the exitflag of the optmization, and
        % the value of the slack variables associated to soft constraints.
        
        % Solve MPC problem using batch approach
        [H,f,A,b,lb,ub] = obj.toQp();
        [du_opt,cost,exitflag] = obj.qpSolver.solve(H,f,A,b,lb,ub);
        if isempty(du_opt); error('Solver failed! No solution.'); end
        du_opt = obj.unScale(du_opt);
        
        % Separate the slack variables from the control inputs
        Ns = obj.getNumberSlackVar();
        Np_MPC = obj.Np;
        slack  = du_opt(end-Ns*Np_MPC+1:end);
        du_opt = du_opt(1:end-Ns*Np_MPC);
        du_opt = reshape(du_opt, [obj.Nu+obj.Nw obj.Nt]);
        
        % Return the control input
        u_seq_opt = obj.uPrev + du_opt(1:obj.Nu,:);
        u_opt     = u_seq_opt(1:obj.Nu,1);
        
        % Save the current control inputs for next timestep
        obj.uPrev = u_opt;
    end
    
    function reset(obj)
        % Reset the MPC controller (initialize/reset its states)
        obj.uPrev = zeros(obj.Nu,1);
        obj.qpSolver.reset();
    end
    
    function release(obj)
        % Release resources
        obj.qpSolver.release();
    end
end

methods (Access = private)
    function [HBatch,fBatch,ABatch,bBatch,lbBatch,ubBatch] = toQp(obj)
        % Define linearization point as the current states and the 
        % previous control inputs
        uop = [obj.uPrev; obj.wMeas];

        if obj.plant.Ts == 0
            % Discretie the model with c2d
            sysc = ss(obj.plant.A,obj.plant.B,[],[]);
            sysd = c2d(sysc,obj.Ts,'zoh');
            A = sysd.A;
            B = sysd.B;
            C = obj.plant.C;
            D = obj.plant.D;
        else
            A = obj.plant.A;
            B = obj.plant.B;
            C = obj.plant.C;
            D = obj.plant.D;
        end
        
        % Formulate the stage cost as
        %    J_k = [   e_k   ]^T [Q   T] [   e_k   ] + [fx]^T [   e_k   ]
        %          [Delta u_k]   [T^T R] [Delta u_k]   [fu]   [Delta u_k]
        % where e_k = y_k-y_ref is the reference tracking error and
        % Delta u_k is the perturbed control input.
        Qe = obj.cost.Qe;
        Rd = obj.cost.Rd;
        Ru = obj.cost.Ru;
        Q  = C' * Qe * C;
        R  = D' * Qe * D + Rd +  Ru;
        T  = C' * Qe * D;
        dyRef = obj.yRef - obj.plant.yop;
        fx = (-2 * dyRef' * Qe * C)';
        fu = (-2 * dyRef' * Qe * D)' + (2 * uop' * Ru)';
        
        % Constraints on perturbed variable are
        %   Ax dx + Au du + As slack <= bineq - bop
        Ax = obj.constraints.Ax;
        Au = obj.constraints.Au;
        b  = obj.constraints.b - obj.constraints.bop;
        Ns = obj.getNumberSlackVar();
        As = zeros(size(Ax,1),Ns);
        for i = 1:Ns
            As(obj.softConstraints.index{i},i) = -1;
        end
        
        % Bounds on the perturbed control inputs.
        % Maximum bounds on the known inputs are set to zero to keep
        % the known inputs contant over the prediction horizon
        dumin = [obj.bounds.lb - uop(1:obj.Nu); zeros(obj.Nw,1)];
        dumax = [obj.bounds.ub - uop(1:obj.Nu); zeros(obj.Nw,1)];

        % Get the number of states, inputs, constraints and horizon
        Nx_MPC = size(A,1);         % Number of states
        Nu_MPC = size(B,2);         % Number of control inputs
        Nt_MPC = obj.Nt;            % Control horizon
        Np_MPC = obj.Np;            % Prediction horizon
        Nineq  = size(Ax,1);        % Inequality constraints

        % +-------------------------+
        % |      Cost function      |
        % +-------------------------+
        % Solve the problem by using a batch approch:
        % Write the state equation over the prediction horizon:
        %       x = Sx x0 + Su u
        % where:
        %       x = [x0 ... xN]^T
        %       u = [u0 ... uN-1]^T
        % Since x0 = 0, we have x = Su u

        Su = zeros(Nx_MPC*Np_MPC,Nu_MPC*Np_MPC);
        Su(Nx_MPC+(1:Nx_MPC), 1:Nu_MPC) = B;
        % Fill first column
        for i = 2:Np_MPC-1
            Su(Nx_MPC*i+(1:Nx_MPC), 1:Nu_MPC) = ...
                A * Su(Nx_MPC*(i-1)+(1:Nx_MPC), 1:Nu_MPC);
        end
        % Copy first column to other columns
        for j = 1:Np_MPC-1
            Su(Nx_MPC*(j+1)+1:end, Nu_MPC*j+(1:Nu_MPC)) = ...
                Su(Nx_MPC+1:Nx_MPC*(Np_MPC-j), 1:Nu_MPC);
        end
        % Simplify by using control horizon (this changes the size of
        % Su to a Nx_MPC*Np_MPC x Nu_MPC*Nt_MPC matrix.
        u2uc = eye(Nu_MPC*Np_MPC, Nu_MPC*Nt_MPC);
        for k = Nt_MPC:Np_MPC-1
            u2uc(...
                k*Nu_MPC+(1:Nu_MPC), ...
                (Nt_MPC-1)*Nu_MPC+(1:Nu_MPC)...
            ) = eye(Nu_MPC);
        end
        Su = Su * u2uc;

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

        HBatch = 2 * (Su' * QQ * Su + TT' * Su + Su' * TT + RR);
        fBatch = Su'*Fx + Fu;

        % +-------------------------+
        % |         Bounds          |
        % +-------------------------+
        % Check min <= max
        if sum(dumin > dumax)
            error('Unfeasible bounds on u')
        end

        % Extend matrices dumin, dumax to obtain dulb and duub where:
        %       dulb <= du <= duub
        % with:
        %       du = [du0; du1; ... duN-1]
        ubBatch = zeros(Nu_MPC*Nt_MPC,1);
        for k = 1:Nt_MPC
            ubBatch((k-1)*Nu_MPC+1:k*Nu_MPC,:) = dumax;
        end
        
        lbBatch = zeros(Nu_MPC*Nt_MPC,1);
        for k = 1:Nt_MPC
            lbBatch((k-1)*Nu_MPC+1:k*Nu_MPC,:) = dumin;
        end

        % +-------------------------+
        % | Inequality  constraints |
        % +-------------------------+
        % Inequalities are written as:
        %       Ax dx + Au du <= bineq
        if isempty(Ax) && isempty(Au)
            ABatch = zeros(0,Nu_MPC*Nt_MPC);
            bBatch = [];
        else
            % Create matrix blockdiag(Au,...,Au)
            AAu = zeros(Nineq*Np_MPC, Nu_MPC*Np_MPC);
            for k = 1:Np_MPC
                AAu((k-1)*Nineq+1:k*Nineq,...
                    (k-1)*Nu_MPC+1:k*Nu_MPC) = Au;
            end
            % Use control horizon (this changes the size of
            % AAu to a Nineq*Np_MPC x Nu_MPC*Nt_MPC matrix.
            AAu = AAu * u2uc;

            % Create matrix blockdiag(Ax,...,Ax)
            AAx = zeros(Nineq*Np_MPC, Nx_MPC*Np_MPC);
            for k = 1:Np_MPC
                AAx((k-1)*Nineq+1:k*Nineq,...
                    (k-1)*Nx_MPC+1:k*Nx_MPC) = Ax;
            end

            % Write inequalities as: A_ineq_batch * u <= b_ineq_batch
            ABatch = AAx * Su + AAu;
            bBatch = zeros(Nineq*Np_MPC,1);
            for k = 1:Np_MPC
                bBatch((k-1)*Nineq+1:k*Nineq) = ...
                    b;
            end
        end

        % +-----------------------------------+
        % |    Soft inequality constraints    |
        % +-----------------------------------+
        % Compute the augmented matrix 
        %   [ As  ]
        %   [ ... ]
        %   [ As  ]
        % for inequality constraints
        if Ns > 0
            AAs = zeros(size(ABatch,1), Ns*Np_MPC);
            for k = 1:Np_MPC
                AAs((k-1)*size(As,1)+1:k*size(As,1),...
                    (k-1)*Ns        +1:k*Ns) = ...
                    As;
            end
            
            Hs = zeros(Ns*Np_MPC);
            for k = 1:Np_MPC
                Hs((k-1)*Ns+(1:Ns), (k-1)*Ns+(1:Ns)) = ...
                    diag(obj.softConstraints.weight);
            end

            % Modify the optimization problem
            HBatch = blkdiag(HBatch, Hs);
            fBatch = [fBatch; zeros(Ns*Np_MPC,1)];
            ABatch = [ABatch AAs];
            lbBatch = [lbBatch;  zeros(Ns*Np_MPC,1)];
            ubBatch = [ubBatch;  inf(Ns*Np_MPC,1)];
        end

        % +-------------------------+
        % |   Scale the QP problem  |
        % +-------------------------+
        if ~isempty(obj.scaling.input) && ...
                ~(isempty(obj.scaling.slack) && (Ns > 0))
            Ddiag = ones(size(HBatch,1),1);
            for k = 1:Nt_MPC
                Ddiag((k-1)*Nu_MPC+(1:Nu_MPC)) = obj.scaling.input;
            end
            if Ns > 0
                for k = 1:Np_MPC
                    Ddiag(end-Ns*k+1:end-Ns*(k-1)) = obj.scaling.slack;
                end
            end
            %Ddiag = 1./Ddiag;
            obj.scaling.factors = Ddiag;
            D = diag(Ddiag);
            HBatch  = D' * HBatch * D;
            fBatch  = Ddiag .* fBatch;
            ABatch  = ABatch * D;
            lbBatch = Ddiag .\ lbBatch;
            ubBatch = Ddiag .\ ubBatch;
        end
        
        % Make sure H is symmetric
        HBatch = (HBatch+HBatch') / 2;
    end

    function du = unScale(obj,du)
        if ~isempty(obj.scaling.factors)
            Ddiag = obj.scaling.factors;
            du = Ddiag .* du;
        end
    end
    
    function Ns = getNumberSlackVar(obj)
        % Return the number of slack variables
        if isempty(obj.softConstraints)
            Ns = 0;
        else
            Ns = numel(obj.softConstraints.index);
        end
    end
end

end

