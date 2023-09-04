function FcnOut=Fcn4ImplicitEuler(xIn)
%
global xi;
h=0.001;
alpha = 25; beta = 50; %经验值5~50
%alpha = 0; beta = 0; %
q=xIn(1:3);v=xIn(4:6);
%
A = [1 0 0; 0 1 0; 0 0 0];
B = [0; -9.8;0];
P(1,1) = -v(3) ^ 2 * cos(q(3)); P(2,1) = v(3) ^ 2 * sin(q(3)); %加速度约束
phi = [q(1) - sin(q(3)); q(2) + cos(q(3))];
phiq = [1 0 -cos(q(3)); 0 1 -sin(q(3))];
phiT = phiq * v;
P1 = P - 2 * alpha * phiT - (beta ^ 2) * phi;
LEFT = [A phiq'; phiq zeros(2)];
RIGHT = [B; P1];
X = LEFT\RIGHT;
a=X(1:3);
delta=[v;a];

FcnOut=[xi+h*delta-[q;v];
       phiT];
end