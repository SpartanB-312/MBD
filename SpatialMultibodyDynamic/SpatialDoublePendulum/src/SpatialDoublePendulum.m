%% Spatial Double Pendulum
%
%% 初值及测试
% A = [1 0 0; 0 1 0; 0 0 1];
% B = [0; -9.8; -9.8 * sin(pi / 6)];
m1 = 1; m2 = 1;
M = diag([m1, m1, m1, m2, m2, m2]);
L1 = 2; r1 = 0.1;
L2 = 2; r2 = 0.1;
J1 = diag([1/12 * m1 * L1 ^ 2; 1/12 * m1 * L1 ^ 2; 1/2 * m1 * r1 ^ 2]);
J2 = diag([1/12 * m2 * L2 ^ 2; 1/12 * m2 * L2 ^ 2; 1/2 * m2 * r2 ^ 2]);
J = diag([1/12 * m1 * L1 ^ 2, 1/12 * m1 * L1 ^ 2, 1/2 * m1 * r1 ^ 2, 1/12 * m2 * L2 ^ 2, 1/12 * m2 * L2 ^ 2, 1/2 * m2 * r2 ^ 2]);

F1 = m1 * [0; 0; -9.8];
F2 = m2 * [0; 0; -9.8];
F = [F1; F2];
N1 = [0; 0; 0];
N2 = [0; 0; 0];
N = [N1; N2];

q10 = [0; 0.5 * L1 / 2; -sqrt(3) / 2 * L1 / 2; pi / 6; 0; 0]; v10 = [0; 0; 0; 0; 0; 0]; %初始状态
R1 = CgaCal(q10(4, 1), q10(5, 1), q10(6, 1));

q20 = [q10(1:3) + R1 * [0; 0; -L2 / 2] + [0; 0.5 * L2 / 2; -sqrt(3) / 2 * L2 / 2]; pi / 6; 0; 0]; v20 = [0; 0; 0; 0; 0; 0]; %初始状态
R2 = CgaCal(q20(4, 1), q20(5, 1), q20(6, 1));

Phi = [q10(1:3) + R1 * [0; 0; L1 / 2];
                        -q10(1:3) - R1 * [0; 0; -L2 / 2] + q20(1:3) + R2 * [0; 0; L2 / 2]]; %约束
%Phipi = dCgaCal(q0(4,1),q0(5,1),q0(6,1),[0;0;1]);
Phipi = [-R1 * Vec2Mat([0; 0; L1 / 2]), zeros(3);
                        R1 * Vec2Mat([0; 0; -L1 / 2]), -R2 * Vec2Mat([0; 0; L2 / 2])];
Phixyz = [1 0 0 0 0 0;
          0 1 0 0 0 0;
          0 0 1 0 0 0;
          -1 0 0 1 0 0
          0 -1 0 0 1 0;
          0 0 -1 0 0 1];
%Phiq=[Phixyz,Phipi];%雅可比
h = 0.001; T = 5; t = 0:h:T; %步长及时间
% global xi;
% xi=[q0;v0];
column = zeros(6, (length(t) - 1));
q1 = q10; v1 = v10;
a1 = zeros(6, length(t));
q2 = q20; v2 = v20;
a2 = zeros(6, length(t));
E = eye(6);
alpha = 25; beta = 50; %经验值5~50
%alpha=1/h;beta=sqrt(2)/h;
w1_ = [v1(4, 1); v1(5, 1); v1(6, 1)];
wMat1_ = Vec2Mat(w1_);
w2_ = [v2(4, 1); v2(5, 1); v2(6, 1)];
wMat2_ = Vec2Mat(w2_);
P = [-R1 * wMat1_ * wMat1_ * [0; 0; L1 / 2];
                              R1 * wMat1_ * wMat1_ * [0; 0; -L1 / 2] - R2 * wMat2_ * wMat2_ * [0; 0; L2 / 2]];
PhiT = Phixyz * [v1(1:3, 1); v2(1:3, 1)] + Phipi * [v1(4:6, 1); v2(4:6, 1)]; %对时间全微分
P1 = P - 2 * alpha * PhiT - beta ^ 2 * Phi; %稳定形式
LHS = [M, zeros(6), Phixyz';
       zeros(6), J, Phipi';
       Phixyz, Phipi, zeros(6)]; %左端系数
RHS = [F1; F2; R1 * N1 - wMat1_ * J1 * w1_; R2 * N2 - wMat2_ * J2 * w2_; P1]; %右端
X = (LHS ^ -1) * RHS; %上三行加速度；下两行lamda
%% 组装
q1 = [q1 column]; v1 = [v1 column];
q2 = [q2 column]; v2 = [v2 column];
%%
for i = 1:(length(t) - 1)
    R1 = CgaCal(q1(4, i), q1(5, i), q1(6, i));
    R2 = CgaCal(q2(4, i), q2(5, i), q2(6, i));
    D1 = AngleRateDCal(q1(4, i), q1(5, i));
    D2 = AngleRateDCal(q2(4, i), q2(5, i));
    Phi = [q1(1:3, i) + R1 * [0; 0; L1 / 2];
           -q1(1:3, i) - R1 * [0; 0; -L2 / 2] + q2(1:3, i) + R2 * [0; 0; L2 / 2]]; %约束
    Phipi = [-R1 * Vec2Mat([0; 0; L1 / 2]), zeros(3);
             R1 * Vec2Mat([0; 0; -L1 / 2]), -R2 * Vec2Mat([0; 0; L2 / 2])];
    Phixyz = [1 0 0 0 0 0;
              0 1 0 0 0 0;
              0 0 1 0 0 0;
              -1 0 0 1 0 0
              0 -1 0 0 1 0;
              0 0 -1 0 0 1];
    w1_ = [v1(4, 1); v1(5, 1); v1(6, 1)];
    wMat1_ = Vec2Mat(w1_);
    w2_ = [v2(4, 1); v2(5, 1); v2(6, 1)];
    wMat2_ = Vec2Mat(w2_);
    P = [-R1 * wMat1_ * wMat1_ * [0; 0; L1 / 2];
         R1 * wMat1_ * wMat1_ * [0; 0; -L1 / 2] - R2 * wMat2_ * wMat2_ * [0; 0; L2 / 2]];
    PhiT = Phixyz * [v1(1:3, i); v2(1:3, i)] + Phipi * [v1(4:6, i); v2(4:6, i)]; %对时间全微分
    P1 = P - 2 * alpha * PhiT - beta ^ 2 * Phi; %稳定形式
    LHS = [M, zeros(6), Phixyz';
           zeros(6), J, Phipi';
           Phixyz, Phipi, zeros(6)]; %左端系数
    RHS = [F1; F2; R1 * N1 - wMat1_ * J1 * w1_; R2 * N2 - wMat2_ * J2 * w2_; P1]; %右端
    X = LHS \ RHS;
    a1(:, i) = [X(1:3); X(7:9)];
    a2(:, i) = [X(4:6); X(10:12)];
    v1(:, i + 1) = v1(:, i) + h * a1(:, i);
    v2(:, i + 1) = v2(:, i) + h * a2(:, i);
    q1(:, i + 1) = q1(:, i) + h * [eye(3),zeros(3);zeros(3),D1] * v1(:, i);
    q2(:, i + 1) = q2(:, i) + h * [eye(3),zeros(3);zeros(3),D2] * v2(:, i);
end

%% Figure
figure(1),
plot(t,q1(4,:)*180/pi,'-r')
hold on,grid on,
plot(t,q2(4,:)*180/pi,'-g')