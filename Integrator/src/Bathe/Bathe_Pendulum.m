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
ddq = zeros(3, length(t));
E = eye(3);
P = zeros(2, 1);
P(1) = -dq(3) ^ 2 * cos(q(3)); P(2) = dq(3) ^ 2 * sin(q(3)); %加速度约束
P(1) = -dq(3) ^ 2 * sin(q(3)); P(2) = dq(3) ^ 2 * cos(q(3)); %加速度约束
phiT = phiq * dq; %对时间全微分
LEFT = [A phiq'; phiq zeros(2)]; %左端系数
RIGHT = [B; P]; %右端
X = (LEFT ^ -1) * RIGHT; %上三行加速度；下两行lamda
ddq(:, 1) = X(1:3);
l=X(4:5);
dq(:, 2) = dq(:, 1) + h* ddq(:, 1);
q(:, 2) = q(:, 1) + h * dq(:, 1);
%% 组装
q = [q column]; dq = [dq column]; l=[l column(1:2,:) zeros(2,1)];
%% Bathe

for i = 1:(length(t) - 1)

    qt=q(:,i)+dq(:,i)*h/2+ddq(:,i)*h^2/16;
    lt=l(:,i);
    detw=ones(1,5);
    j=1;
    %步长中段
    while max(abs(detw))>1e-4
        dqt=4*(qt-q(:,i))/h-dq(:,i);
        ddqt=16*(qt-q(:,i)-h*dq(:,i)/2)/h^2-ddq(:,i);
        phi = [qt(1) - sin(qt(3)); qt(2) + cos(qt(3))];
        phiq = [1 0 -cos(qt(3)); 0 1 -sin(qt(3))];
        phiqlq = [0,0,0;0,0,0;0,0,sin(qt(3))-cos(qt(3))];
%         St = [h^2/4*(A*(16/h^2))+phiqlq phiq'; phiq zeros(2)];  
%         rq=0.25*h^2*(A*ddqt-B)+phiq'*lt;
%         rl=4*phi/h^2;
        St = [(A*(16/h^2))+phiqlq phiq'; phiq zeros(2)];  
        rq=(A*ddqt-B)+phiq'*lt;
        rl=phi;
        RIGHT = [rq; rl];
        detw = -St\RIGHT;

        qt=qt+detw(1:3);
        lt=lt+detw(4:5);

        j=j+1;
        if j>1000
            break;
        end
    end
    dqt=4*(qt-q(:,i))/h-dq(:,i);
    ddqt=16*(qt-q(:,i)-h*dq(:,i)/2)/h^2-ddq(:,i);
    %
    q(:,i+1)=qt+h*dqt/2+h^2*ddqt/16;
    l(:,i+1)=lt;
    j=1;detw=ones(1,5);
    while max(abs(detw))>1e-4
        dq(:,i+1)=q(:,i)/h-4*qt/h+3*q(:,i+1)/h;
        ddq(:,i+1)=dq(:,i)/h-4*dqt/h+3*dq(:,i+1)/h;
        phi = [q(1,i+1) - sin(q(3,i+1)); q(2,i+1) + cos(q(3,i+1))];
        phiq = [1 0 -cos(q(3,i+1)); 0 1 -sin(q(3,i+1))];
        phiqlq = [0,0,0;0,0,0;0,0,sin(q(3,i+1))-cos(q(3,i+1))];
%         St = [h^2/4*(A*(9/h^2))+phiqlq phiq'; phiq zeros(2)];  
%         rq=0.25*h^2*(A*ddq(:,i+1)-B)+phiq'*l(:.i+1);
%         rl=4*phi/h^2;
        St = [(A*(9/h^2))+phiqlq phiq'; phiq zeros(2)];  
        rq=(A*ddq(:,i+1)-B)+phiq'*l(:,i+1);
        rl=phi;
        RIGHT = [rq; rl];
        detw = -St\RIGHT;

        q(:,i+1)=q(:,i+1)+detw(1:3);
        l(:,i+1)=l(:,i+1)+detw(4:5);

        j=j+1;
        if j>1000
            break;
        end
    end
    dq(:,i+1)=q(:,i)/h-4*qt/h+3*q(:,i+1)/h;
    ddq(:,i+1)=dq(:,i)/h-4*dqt/h+3*dq(:,i+1)/h;


end
%% figure
plot(t,q(3, :)*180/pi);