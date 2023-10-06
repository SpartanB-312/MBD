%
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
Alpha = 25; Beta = 50; %经验值5~50
%alpha=1/h;beta=sqrt(2)/h;
P = zeros(2, 1);
P(1) = -v(3) ^ 2 * cos(q(3)); P(2) = v(3) ^ 2 * sin(q(3)); %加速度约束
P(1) = -v(3) ^ 2 * sin(q(3)); P(2) = v(3) ^ 2 * cos(q(3)); %加速度约束
phiT = phiq * v; %对时间全微分
P1 = P - 2 * Alpha * phiT - Beta ^ 2 * phi; %稳定形式
LEFT = [A phiq'; phiq zeros(2)]; %左端系数
RIGHT = [B; P1]; %右端
X = (LEFT ^ -1) * RIGHT; %上三行加速度；下两行lamda
a(1, 1) = X(1); a(2, 1) = X(2); a(3, 1) = X(3);
l=X(4:5);
v(:, 2) = v(:, 1) + h* a(:, 1);
q(:, 2) = q(:, 1) + h * v(:, 1);
%% 组装
q = [q column]; v = [v column];l=[l column(1:2,:) zeros(2,1)];
%%
alpha=0;
gamma=0.5;eta=(0.5+gamma^2)/4;
C=0;
gamma=0.5+C;eta=(1+C^2)/4;
%%

for i = 1:(length(t) - 1)
    a(:,i+1)=a(:,i);
    detw=ones(1,5);
    j=1;
    while max(abs(detw))>1e-4
        q(:, i + 1) = q(:, i) + h * v(:, i) + 0.5* h^2 *((1-2*eta)*a(:,i)+2*eta*a(:,i+1));
        v(:, i + 1) = v(:, i) + h * ((1-gamma)*a(:, i)+gamma*a(:,i+1));
        
        Kt=diag([0,0,l(1,i+1)*sin(q(3,i+1))+l(2,i+1)*sin(q(3, i+1))])*0;
        Ct=0;
        Pt=eta*h^2*Kt;
        
        phi = [q(1, i+1) - sin(q(3, i+1)); q(2, i+1) + cos(q(3, i+1))];
        phiq = [1 0 -cos(q(3, i+1)); 0 1 -sin(q(3, i+1))];
        phiT = phiq * v(:, i+1);
        St = [A+Pt phiq'; phiq zeros(2)];
        rl=phi/eta/h^2;
        rq=A*a(:,i+1)+phiq'*l(:,i+1)-B;%最后一个B是上一个的B，但此处B不变
        RIGHT = [rq; rl];
        
        detw = St\RIGHT;
        a(:,i+1)=a(:,i+1)-detw(1:3);
        l(:,i+1)=l(:,i+1)-detw(4:5);
        

        j=j+1;
        if j>10000
            break;
        end
    end


end

%%

%% Figure
% lengthoft=length(t);
plot(t,q(3, :)*180/pi);
% plot(t,MEA_ANGLE_1Q-q(3, :)'*180/pi);