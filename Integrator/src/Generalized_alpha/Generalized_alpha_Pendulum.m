%% 初值及测试
A = [1 0 0; 0 1 0; 0 0 0];
B = [0; -9.8; 0];%-9.8 * sin(pi / 6)];
q0 = [0.5; -sqrt(3) / 2; pi / 6]; v0 = [0; 0; 0]; %初始状态
phi = [q0(1) - sin(pi / 6); q0(2) + cos(pi / 6)]; %约束
phiq = [1 0 -cos(pi / 6); 0 1 -sin(pi / 6)]; %雅可比
h = 0.001; T = 5; t = 0:h:T; %步长及时间

column = zeros(3, (length(t) - 2));
q = q0;
v = v0;
a = zeros(3, length(t));
E = eye(3);
alpha = 25; beta = 50; %经验值5~50
%alpha=1/h;beta=sqrt(2)/h;
P = zeros(2, 1);
P(1) = -v(3) ^ 2 * cos(q(3)); P(2) = v(3) ^ 2 * sin(q(3)); %加速度约束
phiT = phiq * v; %对时间全微分
P1 = P - 2 * alpha * phiT - beta ^ 2 * phi; %稳定形式
LEFT = [A phiq'; phiq zeros(2)]; %左端系数
RIGHT = [B; P1]; %右端
X = (LEFT ^ -1) * RIGHT; %上三行加速度；下两行lamda
a(1, 1) = X(1); a(2, 1) = X(2); a(3, 1) = X(3);
v(:, 2) = v(:, 1) + h* a(:, 1);
q(:, 2) = q(:, 1) + h * v(:, 1);
%% 组装
q = [q column]; v = [v column];
%% 广义alpha
%后缀加i与Baumgrate区分
alphai=0;%α在[0,1/3]
gammai=0+alphai;betai=0.25*(0.5+gammai)^2;

for i = 2:(length(t) - 1)
    %B(3, 1) = -9.8 * sin(q(3, i));
    P(1) = -v(3, i) ^ 2 * cos(q(3, i)); P(2) = v(3, i) ^ 2 * sin(q(3, i)); %加速度约束
    phi = [q(1, i) - sin(q(3, i)); q(2, i) + cos(q(3, i))];
    phiq = [1 0 -cos(q(3, i)); 0 1 -sin(q(3, i))];
    phiT = phiq * v(:, i);
    P1 = P - 2 * alpha * phiT - (beta ^ 2) * phi;
    LEFT = [A phiq'; phiq zeros(2)];
    RIGHT = [B; P1];
    X = (LEFT ^ -1) * RIGHT;
    a(1, i) = X(1); a(2, i) = X(2); a(3, i) = X(3);
    v(:, i + 1) = v(:, i) + h*(1-gammai) * a(:, i)+h*gammai*a(:,i-1);
    q(:, i + 1) = q(:, i) + h * v(:, i)+(0.5-betai)*(h^2)*a(:,i-1)+betai*(h^2)*a(:,i-1);

end