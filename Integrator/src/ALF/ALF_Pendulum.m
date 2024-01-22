%% 初值及测试
M = [1 0 0; 0 1 0; 0 0 0];
Q = [0; -9.8; 0];%-9.8 * sin(pi / 6)];
q0 = [0.5; -sqrt(3) / 2; pi / 6]; dq0 = [0; 0; 0]; %初始状态
Phi = [q0(1) - sin(pi / 6); q0(2) + cos(pi / 6)]; %约束
Phiq = [1 0 -cos(pi / 6); 0 1 -sin(pi / 6)]; %雅可比
h = 0.001; T = 5; t = 0:h:T; %步长及时间

column = zeros(3, (length(t) - 2));
q = q0;
dq = dq0;
ddq = zeros(3, length(t));
E = eye(3);
P = zeros(2, 1);
P(1) = -dq(3) ^ 2 * cos(q(3)); P(2) = dq(3) ^ 2 * sin(q(3)); %加速度约束
P(1) = -dq(3) ^ 2 * sin(q(3)); P(2) = dq(3) ^ 2 * cos(q(3)); %加速度约束
PhiT = Phiq * dq; %对时间全微分
LEFT = [M Phiq'; Phiq zeros(2)]; %左端系数
RIGHT = [Q; P]; %右端
X = (LEFT ^ -1) * RIGHT; %上三行加速度；下两行lamda
ddq(:, 1) = X(1:3);
l=X(4:5);
dq(:, 2) = dq(:, 1) + h* ddq(:, 1);
q(:, 2) = q(:, 1) + h * dq(:, 1);
%% 组装
q = [q column]; dq = [dq column]; l=[l column(1:2,:) zeros(2,1)];
%% ALF
alphak=10^6;
alpha=diag(repmat(alphak,1,2));
omegak=10;
Omega=diag(repmat(omegak,1,2));
beta=alpha*Omega^2;
uk=1;
u=diag(repmat(uk,1,2));
for i = 1:(length(t) - 1)

    q(:,i+1)=q(:,i)+h*dq(:,i)+1/2*h^2*ddq(:,i);
    dqb=dq(:,i)+2/h*q(:,i);
    ddqb=ddq(:,i)+4/h*dq(:,i)+4/h^2*q(:,i);
    dq(:,i+1)=2/h*q(:,i+1)-dqb;
    ddq(:,i+1)=4/h^2*q(:,i+1)-ddqb;

    
    Qq=zeros(3,1);
    Phi = [q(1,i+1) - sin(q(3,i+1)); q(2,i+1) + cos(q(3,i+1))];
    Phiq = [1 0 -cos(q(3,i+1)); 0 1 -sin(q(3,i+1))];
    dPhi=Phiq*dq(:,i+1)+zeros(2,1);
    P(1) = -dq(3,i+1) ^ 2 * cos(q(3,i+1)); P(2) = dq(3,i+1) ^ 2 * sin(q(3,i+1));
    P(1) = -dq(3,i+1) ^ 2 * sin(q(3,i+1)); P(2) = dq(3,i+1) ^ 2 * cos(q(3,i+1));
    ddPhi=Phiq*ddq(:,i+1)-P;
%     Phiqqp(:,:,1)=zeros(2,3)';
%     Phiqqp(:,:,2)=zeros(2,3)';
%     Phiqqp(:,:,3)=[0 0 sin(q(3,i+1));0 0 -cos(q(3,i+1))]';
%     temp1=beta*Phi;
    X1=[1 0 -cos(q(3,i+1));0 1 -sin(q(3,i+1));-cos(q(3,i+1)) -sin(q(3,i+1)) q(1,i+1)*sin(q(3,i+1))-q(2,i+1)*cos(q(3,i+1))];
    J=4/h^2*M-Qq+X1.*beta(1)+Phiq'*beta*Phiq;


    sigma=beta*dPhi;
    k=beta*ddPhi;

    detw=ones(1,3);
    j=1;

    while max(abs(detw))>1e-4
%         dqb=dq(:,i)+2/h*q(:,i);
%         ddqb=ddq(:,i)+4/h*dq(:,i)+4/h^2*q(:,i);
%         dq(:,i+1)=2/h*q(:,i+1)-dqb;
%         ddq(:,i+1)=4/h^2*q(:,i+1)-ddqb;

        qk=q(:,i+1);dqk=dq(:,i+1);ddqk=ddq(:,i+1);

        Phi = [q(1,i+1) - sin(q(3,i+1)); q(2,i+1) + cos(q(3,i+1))];
        Phiq = [1 0 -cos(q(3,i+1)); 0 1 -sin(q(3,i+1))];
%         Phiqq(:,:,1)=zeros(2,3);
%         Phiqq(:,:,2)=zeros(2,3);
%         Phiqq(:,:,3)=[0 0 sin(q(3,i+1));0 0 -cos(q(3,i+1))];
        dPhi=Phiq*dq(:,i+1)+zeros(2,1);
        P(1) = -dq(3,i+1) ^ 2 * cos(q(3,i+1)); P(2) = dq(3,i+1) ^ 2 * sin(q(3,i+1));
        P(1) = -dq(3,i+1) ^ 2 * sin(q(3,i+1)); P(2) = dq(3,i+1) ^ 2 * cos(q(3,i+1));
        ddPhi=Phiq*ddq(:,i+1)-P;

        Y=M*ddq(:,i+1)+Phiq'*alpha*(ddPhi+2*Omega*u*dPhi+Omega^2*Phi)-Q;
        Yk=Y;
        detw=-J\Y;
        q(:,i+1)=qk+detw;
        dq(:,i+1)=2/h*q(:,i+1)-dqb;
        ddq(:,i+1)=4/h^2*q(:,i+1)-ddqb;

        dPhi=Phiq*dq(:,i+1)+zeros(2,1);
%         P(1) = -dq(3,i+1) ^ 2 * cos(q(3,i+1)); P(2) = dq(3,i+1) ^ 2 * sin(q(3,i+1));
%         ddPhi=Phiq*ddq(:,i+1)-P;
        
%         sigma=sigma+beta*dPhi;
%         k=k+beta*ddPhi;
% 
%         dq(:,i+1)=(M+Phiq'*beta*Phiq)\(M*dq(:,i+1)-Phiq'*sigma);
%         P(1) = -dq(3,i+1) ^ 2 * cos(q(3,i+1)); P(2) = dq(3,i+1) ^ 2 * sin(q(3,i+1));
%         P(1) = -dq(3,i+1) ^ 2 * sin(q(3,i+1)); P(2) = dq(3,i+1) ^ 2 * cos(q(3,i+1));
%         ddPhi=Phiq*ddq(:,i+1)-P;
%         ddq(:,i+1)=(M+Phiq'*beta*Phiq)\(M*ddq(:,i+1)-Phiq'*(-beta*P+k));
        
        sigma=sigma+beta*dPhi;
        k=k+beta*ddPhi;

        Y=M*ddq(:,i+1)+Phiq'*alpha*(ddPhi+2*Omega*u*dPhi+Omega^2*Phi)-Q;

        yk=Y-Yk;
        sk=q(:,i+1)-qk;
        J=J+((yk-J*sk)*sk')/(sk'*sk);

        j=j+1;
        if j>1000
            break;
        end
    end


end
%% figure
plot(t,q(3, :)*180/pi);