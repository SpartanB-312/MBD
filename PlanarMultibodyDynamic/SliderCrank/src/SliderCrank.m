%
%% 参数
%
g = -9.8;
m1 = 1; l1 = 1; J1 = m1 * (l1 ^ 2)/3;
m2 = 1; l2 = 2; J2 = m2 * (l2 ^ 2)/3;
m3 = 1;
w=pi/3*100;
k=10;c=5;
f=0;
%%
A = diag([m1, m1, J1, m2, m2, J2, m3, m3]);
B = [0; m1 * g; 0; 0; m2 * g; 0; f; m3 * g];
q0 = [l1; 0; 0; 2 * l1 + l2; 0; 0; 2 * l1 + 2 * l2; 0];
%q0 = [-l1; 0; pi; -2 * l1 + l2; 0; 0; -2 * l1 + 2 * l2; 0];
%v0 = [0; 0; 0; 0; 0; 0; 0; 0]; %初始状态
v0 = [-l1*sin(q0(3))*w;
          l1*cos(q0(3))*w;
          w;
          -2*l1*sin(q0(3))*w-l2*sin(q0(6))*2*l1*cos(q0(3))*w/(2*l2*cos(q0(6)));
          l1*cos(q0(3))*w;
          -l1*cos(q0(3))*w/(l2*cos(q0(6))); 
          -2*l1*sin(q0(3))*w-l2*sin(q0(6))*2*l1*cos(q0(3))*w/(2*l2*cos(q0(6)))+l2*sin(q0(6))*2*l1*cos(q0(3))*w/(2*l2*cos(q0(6)));
          0];
q = q0; v = v0;
phi = [q0(1) - l1 * cos(q0(3));
       q0(2) - l1 * sin(q0(3));
       -2 * l1 * cos(q0(3)) + q0(4) - l2 * cos(q0(6));
       -2 * l1 * sin(q0(3)) + q0(5) - l2 * sin(q0(6));
       q0(5)+l2*sin(q0(6));
       %-2 * l1 * sin(q0(3)) + 2 * l2 * sin(q0(6));
       -l2 * cos(q0(6)) - q0(4) + q0(7);
       q0(8)
       q0(3)-w*0]; %约束
phiq = [1 0 l1 * sin(q0(3)) 0 0 0 0 0;
        0 1 -l1 * cos(q0(3)) 0 0 0 0 0;
        0 0 2 * l1 * sin(q0(3)) 1 0 l2 * sin(q0(6)) 0 0;
        0 0 -2 * l1 * cos(q0(3)) 0 1 -l2 * cos(q0(6)) 0 0;
        %0 0 2 * l1 * cos(q0(3)) 0 0 2 * l2 * cos(q0(6)) 0 0;
        0 0 0 -1 0 l2 * sin(q0(6)) 1 0;
        0 0 0 0 1 l2*cos(q0(6)) 0 0;
        0 0 0 0 0 0 0 1;
        0 0 1 0 0 0 0 0]; %雅可比
P = [-l1 * cos(q0(3)) * v0(3) ^ 2;
   -l1 * sin(q0(3)) * v0(3) ^ 2;
   -l1 * cos(q0(3)) * v0(3) ^ 2 - l2 * cos(q0(6)) * v0(6) ^ 2;
   -l1 * sin(q0(3)) * v0(3) ^ 2 - l2 * sin(q0(6)) * v0(6) ^ 2;
   l2*sin(q0(6))*v0(6)^2
   %2 * l1 * sin(q0(3)) * v0(3) ^ 2 + 2 * l2 * sin(q0(6)) * v0(6) ^ 2;
   l2 * cos(q0(6)) * v0(6) ^ 2;
   0;
   0];
h = 0.0001; T = 0.1; t = 0:h:T; %步长及时间

phiT = phiq * v; %对时间全微分
alpha = 25; beta = 50; %经验值5~50
P1 = P - 2 * alpha * phiT - beta ^ 2 * phi; %稳定形式
LEFT = [A phiq'; phiq zeros(length(P))]; %左端系数
RIGHT = [B; P1]; %右端
X = LEFT\RIGHT;

%% 组装
q = q0; v = v0;
column = zeros(8, (length(t) - 1));
q = [q column]; v = [v column];
a = zeros(8, length(t));
%%
for i = 1:length(t) - 1
    phi = [q(1, i) - l1 * cos(q(3, i));
           q(2, i) - l1 * sin(q(3, i));
           -2 * l1 * cos(q(3, i)) + q(4, i) - l2 * cos(q(6, i));
           -2 * l1 * sin(q(3, i)) + q(5, i) - l2 * sin(q(6, i));
           %-2 * l1 * sin(q(3, i)) + 2 * l2 * sin(q(6, i))
           l2 * cos(q(6, i)) + q(4, i) - q(7, i);
           q(5,i)+l2*sin(q(6,i));
           q(8, i);
           q(3,i)-w*h*(i-1)]; %约束
    phiq = [1 0 l1 * sin(q(3, i)) 0 0 0 0 0;
            0 1 -l1 * cos(q(3, i)) 0 0 0 0 0;
            0 0 2 * l1 * sin(q(3, i)) 1 0 l2 * sin(q(6, i)) 0 0;
            0 0 -2 * l1 * cos(q(3, i)) 0 1 -l2 * cos(q(6, i)) 0 0;
            %0 0 2 * l1 * cos(q(3, i)) 0 0 2 * l2 * cos(q(6, i)) 0 0;
            0 0 0 1 0 -l2 * sin(q(6, i)) -1 0;
            0 0 0 0 1 l2*cos(q(6,i)) 0 0;
            0 0 0 0 0 0 0 1
            0 0 1 0 0 0 0 0]; %雅可比
    P = [-l1 * cos(q(3, i)) * v(3, i) ^ 2;
       -l1 * sin(q(3, i)) * v(3, i) ^ 2;
       -2*l1 * cos(q(3, i)) * v(3, i) ^ 2 - l2 * cos(q(6, i)) * v(6, i) ^ 2;
       -2*l1 * sin(q(3, i)) * v(3, i) ^ 2 - l2 * sin(q(6, i)) * v(6, i) ^ 2;
        %l1 * sin(q(3, i)) * v(3, i) ^ 2 + 2 * l2 * sin(q(6, i)) * v(6, i) ^ 2;
       l2 * cos(q(6, i)) * v(6, i) ^ 2;
       l2*sin(q(6,i))*v(6,i)^2;
       0;
       0];
    phit=[0;0;0;0;0;0;0;-w];
    phiT = phiq * v(:, i)+phit;
    P1 = P - 2 * alpha * phiT - (beta ^ 2) * phi;
    LEFT = [A phiq'; phiq zeros(length(P))];
    RIGHT = [B; P1];
    X = LEFT \ RIGHT;
    a(:, i) = X(1:8);
    v(:, i + 1) = v(:, i) + h * a(:, i);
    q(:, i + 1) = q(:, i) + h * v(:, i);
end

%% figure
figure(1);
hold on;
grid on;
plot(t, q(7, :),'-g');
plot(t, q(8, :),'-r');
legend('滑块X方向','滑块Y方向');

figure(2);
hold on;
grid on;
plot(t, v(7, :),'-g');
legend('滑块X方向速度');

figure(3);
hold on;
grid on;
plot(t, a(7, :),'-g');
legend('滑块X方向加速度');