% UMVINIT Initialize the model parameters.

A = [0.9550    0.2250  135.7879;
     0.0450    0.7750 -135.7879;
    -0.0000    0.0001    0.9893];
B = [0;
     0;
     0];
H = [0.0012;
     0.0088;
     0.0000];
C = [0.0001         0         0;
     0.0005   -0.0025   -1.5791];
D = [0;
     0];
G = [0;
     1];
Ts = 0.01;
stateInit = [10; 0; 0];
stateInitUmv = [10; 0; 0];
Q = eye(size(A,1))*1e-3;
R = eye(size(C,1))*1e-3;
PInit = eye(size(A,1))*1e-1;