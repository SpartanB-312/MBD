%% 初值及测试
A = [1 0 0; 0 1 0; 0 0 0];
B = [0; -9.8; 0];%-9.8 * sin(pi / 6)];
q0 = [0.5; -sqrt(3) / 2; pi / 6]; dq0 = [0; 0; 0]; %初始状态
phi = [q0(1) - sin(pi / 6); q0(2) + cos(pi / 6)]; %约束
phiq = [1 0 -cos(pi / 6); 0 1 -sin(pi / 6)]; %雅可比
h = 0.001; T = 5; t = 0:h:T; %步长及时间

column = zeros(3, (length(t) - 2));
q = q0;
dq = dq0;
a = zeros(3, length(t));
E = eye(3);
Alpha = 25; Beta = 50; %经验值5~50
%alpha=1/h;beta=sqrt(2)/h;
P = zeros(2, 1);
P(1) = -dq(3) ^ 2 * cos(q(3)); P(2) = dq(3) ^ 2 * sin(q(3)); %加速度约束
P(1) = -dq(3) ^ 2 * sin(q(3)); P(2) = dq(3) ^ 2 * cos(q(3)); %加速度约束
phiT = phiq * dq; %对时间全微分
LEFT = [A phiq'; phiq zeros(2)]; %左端系数
RIGHT = [B; P]; %右端
X = (LEFT ^ -1) * RIGHT; %上三行加速度；下两行lamda
ddq(1, 1) = X(1); ddq(2, 1) = X(2); ddq(3, 1) = X(3);
l=X(4:5);
dq(:, 2) = dq(:, 1) + h* ddq(:, 1);
q(:, 2) = q(:, 1) + h * dq(:, 1);
%% 组装
q = [q column]; dq = [dq column]; l=[l column(1:2,:) zeros(2,1)];
%% 广义alpha
alpha=-1/3;%α在[-1/3,0]
p=0.6;
alpham=(2*p-1)/(1+p);alphaf=p/(1+p);
beta=0.25*(1+alphaf-alpham)^2;
betab=(1-alpham)/(h^2*beta*(1-alphaf));
gamma=0.5+alphaf-alpham;
gammab=gamma/h/beta;

a=ddq(:,1);

for i = 1:(length(t) - 1)
    q(:,i+1)=q(:,i)+h*dq(:,i)+h^2*(0.5-beta)*a;
    dq(:,i+1)=dq(:,i)+h*(1-gamma)*a;

    a=(alphaf*ddq(:,i)-alpham*a)/(1-alpham);
    l(:,i+1)=zeros(2,1);

    q(:,i+1)=q(:,i+1)+h^2*beta*a;
    dq(:,i+1)=dq(:,i+1)+h*gamma*a;
    ddq(:,i+1)=zeros(3,1);

    detw=ones(1,5);
    j=1;

    while max(abs(detw))>1e-4

        phi = [q(1, i+1) - sin(q(3, i+1)); q(2, i+1) + cos(q(3, i+1))];
        phiq = [1 0 -cos(q(3, i+1)); 0 1 -sin(q(3, i+1))];
        phiqlq = [0,0,0;0,0,0;0,0,sin(q(3,i+1))-cos(q(3,i+1))];
        Ct=zeros(3,3);
        Kt=phiqlq;
        St = [A*betab+Ct*betab+Kt phiq'; phiq zeros(2)];
        rl=phi;
        rq=A*ddq(:,i+1)+phiq'*l(:,i+1)-B;
        RIGHT = [rq; rl];
        detw = -St\RIGHT;

        q(:,i+1)=q(:,i+1)+detw(1:3);
        dq(:,i+1)=dq(:,i+1)+gammab*detw(1:3);
        ddq(:,i+1)=ddq(:,i+1)+betab*detw(1:3);
        l(:,i+1)=l(:,i+1)+detw(4:5);


        j=j+1;
        if j>10000
            break;
        end
    end

    a=a+(1-alphaf)/(1-alpham)*ddq(:,i+1);

end

%% figure
plot(t,q(3, :)*180/pi);