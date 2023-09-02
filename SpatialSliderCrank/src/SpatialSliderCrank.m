%% Spatial Slider Crank
%[Haug]CAKD-I p396
%%
m1 = 0.5; m2 = 5;
M = diag([m1, m1, m1, m2, m2, m2]);
J1=diag([0.2,0.2,0.2]);
J2=diag([0.2,0.2,0.2]);
J=diag([0.2,0.2,0.2,0.2,0.2,0.2]);
F1 = m1 * [0; 0; -9.8];
F2 = m2 * [0; 0; -9.8];
F = [F1; F2];
N1 = [0; 0; 0];
N2 = [0; 0; 0];
N = [N1; N2];
omega0=-2*pi;

q10 = [0; 0.1; 0.12; pi / 2; 0; 0]; v10 = [0; 0; 0; omega0; 0; 0]; %初始状态
R1 = CgaCal(q10(4, 1), q10(5, 1), q10(6, 1));

%[rll,pth,yaw]=EulerAngleCal([0.15;0;0],[0.1;-0.05;-0.1]);
%q20 = [0.1; 0.05; 0.1; rll; pth; yaw]; v20 = [0; 0; 0; 0; 0; 0]; %初始状态
%R2 = CgaCal(q20(4, 1), q20(5, 1), q20(6, 1));

q20 = [0.2; 0; 0; 0; 0; 0]; v20 = [-0.2513; 0; 0; 0; 0; 0]; %初始状态
R2 = CgaCal(q20(4, 1), q20(5, 1), q20(6, 1));
%Constraints
%crank-ground revolute i:ground
%crank-connecting_rod spherical
%connecting_rod-slider revolute-cylindrical
%slider-ground translational p1,p2,d1 i:ground j:
%connectinig_rod-slider distance i:crank
Phi = [q10(1:3,1)-[0; 0.1; 0.12];
       [0;1;0]'*R1*[1;0;0];
       [0;0;1]'*R1*[1;0;0];
       %-q10(1:3,1)-R1*[0;0.08;0]+q20(1:3,1)+R2*[-0.15;0;0];
       %[0;1;0]'*R2'*R3*[1;0;0];
       %[0;1;0]'*R3'*[0.1;0;0];
       %[0;0;1]'*R3'*[0.1;0;0];
       [0;1;0]'*R2*[1;0;0];
       [0;0;1]'*R2*[1;0;0];
       [0;1;0]'*(q20(1:3,1));
       [0;0;1]'*(q20(1:3,1));
       [0;1;0]'*R2*[0;0;1];
       %(q20(1:3,1)-q30(1:3,1))'*(q20(1:3,1)-q30(1:3,1))
       (q20(1:3,1)-q10(1:3,1)-R1*[0;0.08;0])'*(q20(1:3,1)-q10(1:3,1)-R1*[0;0.08;0])-0.09;
       q10(4,1)-q10(4,1)-omega0*0]; %约束

Phipi = [-R1*Vec2Mat([0;0;0]),zeros(3,3);
         -[0;1;0]'*R1*Vec2Mat([1;0;0]),zeros(1,3);
         -[0;0;1]'*R1*Vec2Mat([1;0;0]),zeros(1,3);
         zeros(1,3),-[0;1;0]'*R2*Vec2Mat([1;0;0]);
         zeros(1,3),-[0;0;1]'*R2*Vec2Mat([1;0;0]);
         %zeros(1,3),-[0;1;0]*Vec2Mat([0;0;0])+(-q20(1:3,1))*Vec2Mat([0;1;0]);
         %zeros(1,3),-[0;0;1]*Vec2Mat([0;0;0])+(-q20(1:3,1))*Vec2Mat([0;0;1]);
         zeros(1,3),zeros(1,3);
         zeros(1,3),zeros(1,3);
         zeros(1,3),-[0;1;0]'*R2*Vec2Mat([0;0;1]);
         2*(q20(1:3,1)-q10(1:3,1)-R1*[0;0.08;0])'*R1*Vec2Mat([0;0.08;0]),zeros(1,3);
         1,0,0,zeros(1,3)];
Phixyz = [eye(3),zeros(3,3);
          zeros(1,3),zeros(1,3);
          zeros(1,3),zeros(1,3);
          zeros(1,3),zeros(1,3);
          zeros(1,3),zeros(1,3);
          zeros(1,3),[0;1;0]';
          zeros(1,3),[0;0;1]';
          zeros(1,3),zeros(1,3);
          -2*(q20(1:3,1)-q10(1:3,1)-R1*[0;0.08;0])',2*(q20(1:3,1)-q10(1:3,1)-R1*[0;0.08;0])';
          zeros(1,3),zeros(1,3)];
Phit = [0;0;0;0;0;0;0;0;0;0;0;-omega0];
%Phiq=[Phixyz,Phipi];%雅可比
h = 0.001; T = 1; t = 0:h:T; %步长及时间
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
P = [gammas(eye(3),R1,zeros(3,1),w1_,[0;0.1;0.12],zeros(3,1));
     gammad1([0;1;0],[1;0;0],eye(3),R1,zeros(3,1),w1_);
     gammad1([0;0;1],[1;0;0],eye(3),R1,zeros(3,1),w1_);
     gammad1([0;1;0],[1;0;0],eye(3),R2,zeros(3,1),w2_);
     gammad1([0;0;1],[1;0;0],eye(3),R2,zeros(3,1),w2_);
     gammad2([0;1;0],eye(3),R2,zeros(3,1),w2_,zeros(3,1),v2(1:3,1),zeros(3,1),zeros(3,1),(-q20(1:3,1)));
     gammad2([0;0;1],eye(3),R2,zeros(3,1),w2_,zeros(3,1),v2(1:3,1),zeros(3,1),zeros(3,1),(-q20(1:3,1)));
     gammad1([1;0;0],[0;0;1],eye(3),R2,zeros(3,1),w2_);
     gammass(R1,R2,w1_,w2_,v1(1:3,1),v2(1:3,1),[0;0.08;0],zeros(3,1),(q20(1:3,1)-q10(1:3,1)-R1*[0;0.08;0]));
     gammatheta(eye(3),R1,zeros(3,1),w1_,[1;0;0])];
PhiT = Phixyz * [v1(1:3, 1); v2(1:3, 1)] + Phipi * [v1(4:6, 1); v2(4:6, 1)] + Phit; %对时间全微分
P1 = P - 2 * alpha * PhiT - beta ^ 2 * Phi; %稳定形式
LHS = [M, zeros(6), Phixyz';
       zeros(6), J, Phipi';
       Phixyz, Phipi, zeros(12)]; %左端系数
RHS = [F1; F2; R1 * N1 - wMat1_ * J1 * w1_; R2 * N2 - wMat2_ * J2 * w2_; P1]; %右端
X = LHS \ RHS;
%% 组装
q1 = [q1 column]; v1 = [v1 column];
q2 = [q2 column]; v2 = [v2 column];
%%
for i = 1:(length(t) - 1)
    R1 = CgaCal(q1(4, i), q1(5, i), q1(6, i));
    R2 = CgaCal(q2(4, i), q2(5, i), q2(6, i));
    D1 = AngleRateDCal(q1(4, i), q1(5, i));
    D2 = AngleRateDCal(q2(4, i), q2(5, i));
    Phi = [q1(1:3,i)-[0; 0.1; 0.12];
       [0;1;0]'*R1*[1;0;0];
       [0;0;1]'*R1*[1;0;0];
       %-q10(1:3,1)-R1*[0;0.08;0]+q20(1:3,1)+R2*[-0.15;0;0];
       %[0;1;0]'*R2'*R3*[1;0;0];
       %[0;1;0]'*R3'*[0.1;0;0];
       %[0;0;1]'*R3'*[0.1;0;0];
       [0;1;0]'*R2*[1;0;0];
       [0;0;1]'*R2*[1;0;0];
       [0;1;0]'*(-q2(1:3,i));
       [0;0;1]'*(-q2(1:3,i));
       [0;1;0]'*R2*[0;0;1];
       %(q20(1:3,1)-q30(1:3,1))'*(q20(1:3,1)-q30(1:3,1))
       (q2(1:3,i)-q1(1:3,i)-R1*[0;0.08;0])'*(q2(1:3,i)-q1(1:3,i)-R1*[0;0.08;0])-0.09;
       q1(4,i)-q10(4,1)-omega0*h*(i-1)]; %约束

Phipi = [-R1*Vec2Mat([0;0;0]),zeros(3,3);
         -[0;1;0]'*R1*Vec2Mat([1;0;0]),zeros(1,3);
         -[0;0;1]'*R1*Vec2Mat([1;0;0]),zeros(1,3);
         zeros(1,3),-[0;1;0]'*R2*Vec2Mat([1;0;0]);
         zeros(1,3),-[0;0;1]'*R2*Vec2Mat([1;0;0]);
         %zeros(1,3),-[0;1;0]*Vec2Mat([0;0;0])+(-q20(1:3,1))*Vec2Mat([0;1;0]);
         %zeros(1,3),-[0;0;1]*Vec2Mat([0;0;0])+(-q20(1:3,1))*Vec2Mat([0;0;1]);
         zeros(1,3),zeros(1,3);
         zeros(1,3),zeros(1,3);
         zeros(1,3),-[0;1;0]'*R2*Vec2Mat([0;0;1]);
         2*(q2(1:3,i)-q1(1:3,i)-R1*[0;0.08;0])'*R1*Vec2Mat([0;0.08;0]),zeros(1,3)
         1,0,0,zeros(1,3)];
Phixyz = [eye(3),zeros(3,3);
          zeros(1,3),zeros(1,3);
          zeros(1,3),zeros(1,3);
          zeros(1,3),zeros(1,3);
          zeros(1,3),zeros(1,3);
          zeros(1,3),[0;1;0]';
          zeros(1,3),[0;0;1]';
          zeros(1,3),zeros(1,3);
          -2*(q2(1:3,i)-q1(1:3,i)-R1*[0;0.08;0])',2*(q2(1:3,i)-q1(1:3,i)-R1*[0;0.08;0])';
          zeros(1,3),zeros(1,3)];
Phit = [0;0;0;0;0;0;0;0;0;0;0;-omega0];
w1_ = [v1(4, i); v1(5, i); v1(6, i)];
wMat1_ = Vec2Mat(w1_);
w2_ = [v2(4, i); v2(5, i); v2(6, i)];
wMat2_ = Vec2Mat(w2_);
P = [gammas(eye(3),R1,zeros(3,1),w1_,[0;0.1;0.12],zeros(3,1));
     gammad1([0;1;0],[1;0;0],eye(3),R1,zeros(3,1),w1_);
     gammad1([0;0;1],[1;0;0],eye(3),R1,zeros(3,1),w1_);
     gammad1([0;1;0],[1;0;0],eye(3),R2,zeros(3,1),w2_);
     gammad1([0;0;1],[1;0;0],eye(3),R2,zeros(3,1),w2_);
     gammad2([0;1;0],eye(3),R2,zeros(3,1),w2_,zeros(3,1),v2(1:3,i),zeros(3,1),zeros(3,1),(-q2(1:3,i)));
     gammad2([0;0;1],eye(3),R2,zeros(3,1),w2_,zeros(3,1),v2(1:3,i),zeros(3,1),zeros(3,1),(-q2(1:3,i)));
     gammad1([1;0;0],[0;0;1],eye(3),R2,zeros(3,1),w2_);
     gammass(R1,R2,w1_,w2_,v1(1:3,i),v2(1:3,i),[0;0.08;0],zeros(3,1),(q2(1:3,i)-q1(1:3,i)-R1*[0;0.08;0]))
     gammatheta(eye(3),R1,zeros(3,1),w1_,[1;0;0])];
PhiT = Phixyz * [v1(1:3, i); v2(1:3, i)] + Phipi * [v1(4:6, i); v2(4:6, i)] + Phit; %对时间全微分
P1 = P - 2 * alpha * PhiT - beta ^ 2 * Phi; %稳定形式
LHS = [M, zeros(6), Phixyz';
       zeros(6), J, Phipi';
       Phixyz, Phipi, zeros(12)]; %左端系数
RHS = [F1; F2; R1 * N1 - wMat1_ * J1 * w1_; R2 * N2 - wMat2_ * J2 * w2_; P1]; %右端
X = (LHS ^ -1) * RHS; 
    a1(:, i) = [X(1:3); X(7:9)];
    a2(:, i) = [X(4:6); X(10:12)];
    v1(:, i + 1) = v1(:, i) + h * a1(:, i);
    v2(:, i + 1) = v2(:, i) + h * a2(:, i);
    q1(:, i + 1) = q1(:, i) + h * [eye(3),zeros(3);zeros(3),D1] * v1(:, i);
    q2(:, i + 1) = q2(:, i) + h * [eye(3),zeros(3);zeros(3),D2] * v2(:, i);
end

%% Figure
figure(1),
hold on,grid on,
plot(t,q2(1,:),'-g')