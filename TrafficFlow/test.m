Umin   = -10;
Umax   = 10;
deltaU = 0.5;
U      = Umin:deltaU:Umax;

% initial conditions
F0 = zeros(1, length(U));
conditions = {{('null', 0), 1}, {(0, 1), 1 - X(X>0 * X<1)}, {(1, 'null'), 0}};
for i =1:
F0(U <= 0) = -1;
F0(U > 0)  = +1;

% Rankine-Hugoniot condition
speed1 = 