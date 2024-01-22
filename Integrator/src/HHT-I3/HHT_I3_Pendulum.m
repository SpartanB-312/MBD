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
q = [q column]; v = [v column]; l=[l column(1:2,:) zeros(2,1)];
%% 广义alpha
alpha=-1/6;%α在[-1/3,0]
gamma=0.5-alpha;beta=0.25*(0.5+gamma)^2;
p=0.5;
alpham=(2*p-1)/(1+p);alphaf=p/(1+p);
betai=(1-alpham)/(h^2*beta*(1-alphaf));
gammai=gamma/(h*beta);

at=a(:,1);

for i = 1:(length(t) - 1)

%     q(:,i+1)=q(:,i)+h*v(:,i)+0.5*h^2*(1-2*beta)*at;
%     v(:,i+1)=v(:,i)+h*(1-gamma)*at;
%     l(:,i+1)=0;
%     at=1/(1-alpham)*(alphaf*a(:,i)-alpham*at);
%     q(:,i+1)=q(:,i+1)+h^2*beta*at;
%     v(:,i+1)=v(:,i+1)+h*gamma*at;
%     a(:,i+1)=0;
    a(:,i+1)=a(:,i);
    
    detw=ones(1,5);
    j=1;

    while max(abs(detw))>1e-4
%         q(:,i+1)=q(:,i)+h*v(:,i)+0.5*h^2*((1-2*beta)*a(:,i)+2*beta*a(:,i));
%         v(:,i+1)=v(:,i)+h*((1-gamma)*a(:,i)+gamma*a(:,i+1));
        q(:,i+1)=q(:,i)+h*v(:,i)+0.5*h^2*((1-2*beta)*a(:,i)+2*beta*a(:,i+1));
        v(:,i+1)=v(:,i)+h*((1-gamma)*a(:,i)+gamma*a(:,i+1));
        
        Kt=zeros(3);
        Ct=0;
        Pt=beta*h^2*Kt;
        
        phi = [q(1, i+1) - sin(q(3, i+1)); q(2, i+1) + cos(q(3, i+1))];
        phiq = [1 0 -cos(q(3, i+1)); 0 1 -sin(q(3, i+1))];
        phiT = phiq * v(:, i+1);
        P1 = P - 2 * Alpha * phiT - (Beta ^ 2) * phi;
        St = [A/(1+alpha)+Pt phiq'; phiq zeros(2)];
        rl=phi/beta/h^2;
        rq=A*a(:,i+1)/(1+alpha)+phiq'*l(:,i+1)-B-alpha/(1+alpha)*(phiq'*l(:,i)-B);%最后一个B是上一个的B，但此处B不变
        RIGHT = [rq; rl];
        detw = St\RIGHT;
        detq=detw(1:3);
        detl=detw(4:5);

%         q(:,i+1)=q(:,i+1)+detq;
%         v(:,i+1)=v(:,i+1)+gammai*detq;
%         a(:,i+1)=a(:,i+1)+betai*detq;
%         l(:,i+1)=l(:,i+1)+detl;

        a(:,i+1)=a(:,i+1)-detw(1:3);
        l(:,i+1)=l(:,i+1)-detw(4:5);


        j=j+1;
        if j>10000
            break;
        end
    end

    at=at+(1-alphaf)/(1-alpham)*a(:,i+1);

end
%% figure
plot(t,q(3, :)*180/pi);