%% Baumgrate
%
%% 初值及测试
clear
A = [1 0 0; 0 1 0; 0 0 0];
B = [0; -9.81; 0];
g=-9.81;
q0 = [0.5; -sqrt(3) / 2; pi / 6]; v0 = [0; 0; 0]; %初始状态
h = 0.001; T = 5; t = 0:h:T; %步长及时间
column = zeros(3, (length(t) - 2));
q = q0;
v = v0;
a = zeros(3, length(t));
E = eye(3);
alpha = 1e7; Omega=diag([10,10]); miu = 1;
LEFT = A; %左端系数
RIGHT = B; %右端
%X = LEFT \ RIGHT; %上三行加速度；下两行lamda
a(:, 1) =[-4.2435;-2.45;-4.9];
%% 组装
q = [q column]; v = [v column];
%%
for i = 1:(length(t) - 1)
    Energy(i)=-1*g*q(2,i)+0.5*v(:,i)'*A*v(:,i);
    Phi = [q(1, i) - sin(q(3, i)); q(2, i) + cos(q(3, i))];
    Phiq = [1 0 -cos(q(3, i)); 0 1 -sin(q(3, i))];
    dPhiq = [0 0 v(3,i)*sin(q(3,i)); 0 0 -v(3,i)*cos(q(3,i))];
    dPhi = Phiq * v(:, i);
    LEFT = A+Phiq'*alpha*Phiq;
    penalty = Phiq'*alpha*(dPhiq*v(:,i)+2*Omega*miu*dPhi+Omega^2*Phi);
    RIGHT = B - penalty;
    X = LEFT  \ RIGHT;
    
    a(:, i + 1) = X;
    v(:, i + 1) = v(:, i) + h * a(:, i);
    q(:, i + 1) = q(:, i) + h * v(:, i);
    
end

%% 
figure(1)
plot(t, q(3, :)); title('角度');
figure(2)
load('MEA_ANGLE_1Q.mat');
plot(t,q(3, :) * 180 / pi - MEA_ANGLE_1Q');
figure(3)
plot(t(1:5000),Energy);
