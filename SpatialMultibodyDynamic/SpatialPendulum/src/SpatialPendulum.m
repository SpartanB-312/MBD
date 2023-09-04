%% Spatial Pendulum
%
%% 初值及测试
m=1;
M = [m 0 0; 0 m 0; 0 0 m];
L = 2;r = 0.1;
J = diag([1/12 * m * L^2; 1/12 * m * L^2; 1/2 * m *r^2]);
F = m*[0; 0; -9.8];
N = [0; 0; 0];

q0 = [0; 0.5*L/2; -sqrt(3) / 2 * L/2; pi / 6; 0; 0]; v0 = [0; 0; 0; 0; 0; 0]; %初始状态
R = CgaCal(q0(4,1),q0(5,1),q0(6,1));

Phi = q0(1:3)+R*[0;0;L/2]; %约束
%Phipi = dCgaCal(q0(4,1),q0(5,1),q0(6,1),[0;0;1]);
Phipi=-R*Vec2Mat([0;0;L/2]);
Phixyz = [1 0 0;
        0 1 0;
        0 0 1]; 
%Phiq=[Phixyz,Phipi];%雅可比
h = 0.001; T = 5; t = 0:h:T; %步长及时间
global xi;
xi=[q0;v0];
column = zeros(6, (length(t) - 1));
q = q0;
v = v0;
a = zeros(6, length(t));
E = eye(6);
alpha = 25; beta = 50; %经验值5~50
%alpha=1/h;beta=sqrt(2)/h;
w_=[v(4,1);v(5,1);v(6,1)];
wMat_=Vec2Mat(w_);
P = -R*wMat_*wMat_*[0;0;L/2];
PhiT = Phixyz * v(1:3,1)+Phipi*v(4:6,1); %对时间全微分
P1 = P - 2 * alpha * PhiT - beta ^ 2 * Phi; %稳定形式
LHS = [M, zeros(3), Phixyz';
       zeros(3), J, Phipi';
       Phixyz, Phipi, zeros(3)]; %左端系数
RHS = [F; R*N-wMat_*J*w_; P1]; %右端
X = (LHS ^ -1) * RHS; %上三行加速度；下两行lamda
%% 组装
q = [q column]; v = [v column];
%%
for i = 1:(length(t) - 1)
    R = CgaCal(q(4,i),q(5,i),q(6,i));
    D = AngleRateDCal(q(4,i),q(5,i));
    Phi = q(1:3,i)+R*[0;0;L/2]; %约束
    Phipi=-R*Vec2Mat([0;0;L/2]);
    Phixyz = [1 0 0;
                    0 1 0;
                    0 0 1];
    w_=[v(4,i);v(5,i);v(6,i)];
    wMat_=Vec2Mat(w_);
    P = -R*wMat_*wMat_*[0;0;L/2];
    PhiT = Phixyz * v(1:3,1)+Phipi*v(4:6,1); %对时间全微分
    P1 = P - 2 * alpha * PhiT - beta ^ 2 * Phi; %稳定形式
    LHS = [M, zeros(3), Phixyz';
                zeros(3), J, Phipi';
                Phixyz, Phipi, zeros(3)]; %左端系数
    RHS = [F; N-wMat_*J*w_; P1]; %右端
    X = (LHS ^ -1) * RHS; 
    a(:, i) = X(1:6);
    v(:, i + 1) = v(:, i) + h * a(:, i);
    q(:, i + 1) = q(:, i) + h * [eye(3),zeros(3);zeros(3),D] * v(:, i);
end

%% 隐式欧拉
for i = 1:(length(t) - 1)
    xi=[q(:,i);v(:,i)];
    opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt','Display','off');
    f=fsolve(@Fcn4ImplicitEuler,xi,opts);
    %f=fsolve(@Fcn4ImplicitEuler,xi,optimset('Display','off'));
    q(:,i+1)=f(1:3);
    v(:,i+1)=f(4:6);
    phi(:,i+1) = [q(1, i+1) - sin(q(3, i+1)); q(2, i+1) + cos(q(3, i+1))];
end
%% Figure
