%% 初值
% 单摆 拉格朗日
%
dy = zeros(2, 1);
h = 0.001;
t = 0:h:5; n = length(t);
global u v a y0
u = pi / 6; v = 0; a = -9.8 * sin(u) + 0;

%% fsolve
for i = 2:n
    y0 = [u(i - 1) v(i - 1) a(i - 1)]; %fsolve寻找零点的初始位置
    f = fsolve(@PendulumImplicit, y0);
    u(i) = f(1);
    v(i) = f(2);
    a(i) = f(3);
end

%% 隐式欧拉
for i = 2:n
    a(i-1) = -9.8 * sin(u(i - 1));
    v(i) = v(i - 1) + h * a(i-1);
    u(i) = u(i - 1) + h * v(i-1);

    for j = 1:15
        a(i) = -9.8 * sin(u(i));
        v(i) = v(i - 1) + h * a(i);
        u(i) = u(i - 1) + h * v(i);
    end

end

%% Function
function F = PendulumImplicit(y)
    global y0
    h = 0.001;
    g = 9.8;
    F(1) = y(3) - (-g * sin(y(1))); %加速度
    F(2) = y(2) - (y0(2) + h * y(3)); %速度
    F(3) = y(1) - (y0(1) + h * y(2)); %位置
end
