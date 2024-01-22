%% Baumgrate
%
%% 初值及测试
% A = [1 0 0; 0 1 0; 0 0 1];
% B = [0; -9.8; -9.8 * sin(pi / 6)];
A = [1 0 0; 0 1 0; 0 0 0];
B = [0; -9.8; 0];
q0 = [0.5; -sqrt(3) / 2; pi / 6]; v0 = [0; 0; 0]; %初始状态
phi = [q0(1) - sin(pi / 6); q0(2) + cos(pi / 6)]; %约束
phiq = [1 0 -cos(pi / 6); 0 1 -sin(pi / 6)]; %雅可比
h = 0.001; T = 5; t = 0:h:T; %步长及时间
global xi;
xi=[q0;v0];
column = zeros(3, (length(t) - 1));
q = q0;
v = v0;
a = zeros(3, length(t));
E = eye(3);
alpha = 25; beta = 50; %经验值5~50
%alpha=1/h;beta=sqrt(2)/h;
P = zeros(2, 1);
P(1) = -v(3) ^ 2 * cos(q(3)); P(2) = v(3) ^ 2 * sin(q(3)); %加速度约束
P(1) = -v(3) ^ 2 * sin(q(3)); P(2) = v(3) ^ 2 * cos(q(3)); %加速度约束
phiT = phiq * v; %对时间全微分
P1 = P - 2 * alpha * phiT - beta ^ 2 * phi; %稳定形式
LEFT = [A phiq'; phiq zeros(2)]; %左端系数
RIGHT = [B; P1]; %右端
X = (LEFT ^ -1) * RIGHT; %上三行加速度；下两行lamda
%% 组装
q = [q column]; v = [v column];
%%
for i = 1:(length(t) - 1)
    %B(3, 1) = -9.8 * sin(q(3, i));
    P(1) = -v(3, i) ^ 2 * cos(q(3, i)); P(2) = v(3, i) ^ 2 * sin(q(3, i)); %加速度约束
    P(1) = -v(3,i) ^ 2 * sin(q(3,i)); P(2) = v(3,i) ^ 2 * cos(q(3,i)); %加速度约束
    phi = [q(1, i) - sin(q(3, i)); q(2, i) + cos(q(3, i))];
    phiq = [1 0 -cos(q(3, i)); 0 1 -sin(q(3, i))];
    phiT = phiq * v(:, i);
    P1 = P - 2 * alpha * phiT - (beta ^ 2) * phi;
    LEFT = [A phiq'; phiq zeros(2)];
    RIGHT = [B; P1];
    X = (LEFT ^ -1) * RIGHT;
    a(1, i) = X(1); a(2, i) = X(2); a(3, i) = X(3);
    v(:, i + 1) = v(:, i) + h * a(:, i);
    q(:, i + 1) = q(:, i) + h * v(:, i);
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

%% 改进欧拉
for i = 1:(length(t) - 1)
    %B(3, 1) = -9.8 * sin(q(3, i));
    P(1) = -v(3, i) ^ 2 * cos(q(3, i));
    P(2) = v(3, i) ^ 2 * sin(q(3, i)); %加速度约束
    phi = [q(1, i) - sin(q(3, i));
           q(2, i) + cos(q(3, i))];
    phiq = [1 0 -cos(q(3, i));
            0 1 -sin(q(3, i))];
    phiT = phiq * v(:, i);
    P1 = P - 2 * alpha * phiT - (beta ^ 2) * phi;
    LEFT = [A phiq';
            phiq zeros(2)];
    RIGHT = [B; P1];
    X = (LEFT ^ -1) * RIGHT;
    a(1, i) = X(1); a(2, i) = X(2); a(3, i) = X(3);
    v(:, i + 1) = v(:, i) + h * a(:, i);
    q(:, i + 1) = q(:, i) + h * v(:, i);

    for j = 1:10
        %B(3, 1) = -9.8 * sin(q(3, i + 1));
        P(1, 1) = -v(3, i + 1) ^ 2 * cos(q(3, i + 1));
        P(2, 1) = v(3, i + 1) ^ 2 * sin(q(3, i + 1)); %加速度约束
        phi = [q(1, i + 1) - sin(q(3, i + 1));
               q(2, i + 1) + cos(q(3, i))];
        phiq = [1 0 -cos(q(3, i + 1));
                0 1 -sin(q(3, i + 1))];
        phiT = phiq * v(:, i + 1);
        P1 = P - 2 * alpha * phiT - (beta ^ 2) * phi;
        LEFT = [A phiq';
                phiq zeros(2)];
        RIGHT = [B; P1];
        X = (LEFT ^ -1) * RIGHT;
        a(1, i + 1) = X(1); a(2, i + 1) = X(2); a(3, i + 1) = X(3);
        v(:, i + 1) = v(:, i) + h * (a(:, i + 1) + a(:, i)) / 2;
        q(:, i + 1) = q(:, i) + h * (v(:, i + 1) + v(:, i)) / 2;
    end

end

%% Figure
subplot(2, 2, 1), plot(t, q(1, :), t, q(2, :)); title('位置');
subplot(2, 2, 2), plot(t, v(3, :)); title('角速度');
subplot(2, 2, 3), plot(t, a(3, :)); title('角加速度');
subplot(2, 2, 4), plot(t, q(3, :)); title('角度');
