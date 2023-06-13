%% Baumgrate
%
%% 初值及测试
A=[1 0 0 0 0 0;
      0 1 0 0 0 0;
      0 0 1 0 0 0;
      0 0 0 1 0 0
      0 0 0 0 1 0
      0 0 0 0 0 4];
B=[0;-9.8;0;
      0;-9.8;0];
q0=[0.5;-sqrt(3)/2;pi/6;
        1;-sqrt(3);pi/6];
v0=[0;0;0;0;0;0];%初始状态
phi=[q0(1,1) * q0(1,1) + q0(2,1) * q0(2,1) - 1;
         - q0(3,1) - atan( q0(1,1) / q0(2,1) );
         ( q0(1,1) - q0(4,1) )^2 + ( q0(2,1) - q0(5,1))^2 - 1;
         - q0(6,1) - atan( ( q0(4,1) - q0(1,1) ) / ( q0(5,1) - q0(2,1) ) )];%约束
phi(4)=[- q0(6,1) - atan(  q0(4,1)   /  q0(5,1)  )];
phiq=[2*q0(1,1) 2*q0(2,1) 0 0 0 0;
           -q0(2,1) , q0(1,1) ,-1 , 0 , 0 , 0;
           2*( q0(1,1) - q0(4,1) ) , 2*(q0(2,1)-q0(5,1)) , 0 , -2*(q0(1,1)-q0(4,1)) , -2*(q0(2,1)-q0(5,1)) , 0;
           -( q0(5,1)- q0(2,1) ) , ( q0(4,1) - q0(1,1) ) , 0 , -( q0(5,1) - q0(2,1) ) , ( q0(4,1) - q0(1,1) ) , -1];%笛卡尔
phiq(4,1)=0;phiq(4,2)=0;phiq(4,3)=0;phiq(4,4)=-q0(5,1)/4;phiq(4,5)=q0(4,1)/4;
h=0.001;T=5;t=0:h:T;%步长及时间

column=zeros( 6 , (length(t)-1) );
q=q0;
v=v0;
a=zeros( 6 , length(t) );
E=eye(6);
alpha=25;beta=50;
%alpha=1/h;beta=sqrt(2)/h;%经验值5~50


P=-2*[v(1,1) v(2,1) 0 0;
           v(1,1)-v(3,1) -v(1,1)+v(3,1) v(2,1)-v(4,1) -v(2,1)+v(4,1)]*v;%

phiT=phiq*v;
P1=P-2*alpha*phiT-beta^2*phi;%稳定形式
LEFT=[A phiq';phiq zeros(4)];%左端系数
RIGHT=[B;P1];%右端
X=(LEFT^-1)*RIGHT;%上两行加速度；第三行lamda
%% 组装
q=[q column];v=[v column];
ak1=zeros(4,1);ak2=zeros(4,1);ak3=zeros(4,1);ak4=zeros(4,1);
vk1=zeros(4,1);vk2=zeros(4,1);vk3=zeros(4,1);vk4=zeros(4,1);
vq1=zeros(4,1);qk2=zeros(4,1);qk3=zeros(4,1);qk4=zeros(4,1);
%% 显式欧拉
for i=1:(length(t)-1)
    P=2*[v(1,i) v(2,i) 0 0;v(1,i)-v(3,i) -v(1,i)+v(3,i) v(2,i)-v(4,i) -v(2,i)+v(4,i)]*v(:,i);
    phi=[q(1,i)*q(1,i)+q(2,i)*q(2,i)-1;(q(1,i)-q(3,i))^2+(q(2,i)-q(4,i))^2-1];%约束
    phiq=[2*q(1,i) 2*q(2,i) 0 0;2*(q(1,i)-q(3,i)) 2*(q(2,i)-q(4,i)) -2*(q(1,i)-q(3,i)) -2*(q(2,i)-q(4,i))];%笛卡尔
    phiT=phiq*v(:,i);
    P1=P-2*alpha*phiT-(beta^2)*phi;
    LEFT=[A phiq';phiq zeros(2,2)];
    RIGHT=[B;P1];
    X=(LEFT^-1)*RIGHT;
    a(1,i)=X(1);a(2,i)=X(2);a(3,i)=X(3);a(4,i)=X(4);
    v(:,i+1)=v(:,i)+h*a(:,i);
    q(:,i+1)=q(:,i)+h*v(:,i);
    
end
%% 隐式欧拉
for i=1:(length(t)-1)
    P=2*[v(1,i) v(2,i) 0 0;v(1,i)-v(3,i) -v(1,i)+v(3,i) v(2,i)-v(4,i) -v(2,i)+v(4,i)]*v(:,i);
    phi=[q(1,i)*q(1,i)+q(2,i)*q(2,i)-1;(q(1,i)-q(3,i))^2+(q(2,i)-q(4,i))^2-1];%约束
    phiq=[2*q(1,i) 2*q(2,i) 0 0;2*(q(1,i)-q(3,i)) 2*(q(2,i)-q(4,i)) -2*(q(1,i)-q(3,i)) -2*(q(2,i)-q(4,i))];%笛卡尔
    phiT=phiq*v(:,i);
    P1=P-2*alpha*phiT-(beta^2)*phi;
    LEFT=[A phiq';phiq zeros(2,2)];
    RIGHT=[B;P1];
    X=(LEFT^-1)*RIGHT;
    a(1,i)=X(1);a(2,i)=X(2);a(3,i)=X(3);a(4,i)=X(4);
    v(:,i+1)=v(:,i)+h*a(:,i);
    q(:,i+1)=q(:,i)+h*v(:,i);
    
    for j=1:15
        
        if max(abs(phi)) < 1e-5 && max(abs(phiT)) < 1e-5
            break
        end
        
        P=2*[v(1,i+1) v(2,i+1) 0 0;v(1,i+1)-v(3,i+1) -v(1,i+1)+v(3,i+1) v(2,i+1)-v(4,i+1) -v(2,i+1)+v(4,i+1)]*v(:,i+1);
        phi=[q(1,i+1)*q(1,i+1)+q(2,i+1)*q(2,i+1)-1;(q(1,i+1)-q(3,i+1))^2+(q(2,i+1)-q(4,i+1))^2-1];%约束
        phiq=[2*q(1,i+1) 2*q(2,i+1) 0 0;2*(q(1,i+1)-q(3,i+1)) 2*(q(2,i+1)-q(4,i+1)) -2*(q(1,i+1)-q(3,i+1)) -2*(q(2,i+1)-q(4,i+1))];%笛卡尔
        phiT=phiq*v(:,i+1);
        
        P1=P-2*alpha*phiT-(beta^2)*phi;
        LEFT=[A phiq';phiq zeros(2,2)];
        RIGHT=[B;P1];
        X=(LEFT^-1)*RIGHT;
        a(1,i+1)=X(1);a(2,i+1)=X(2);a(3,i+1)=X(3);a(4,i+1)=X(4);
        v(:,i+1)=v(:,i)+h*a(:,i+1);
        q(:,i+1)=q(:,i)+h*v(:,i+1);
        
    end
    
%     P=2*[v(1,i+1) v(2,i+1) 0 0;v(1,i+1)-v(3,i+1) -v(1,i+1)+v(3,i+1) v(2,i+1)-v(4,i+1) -v(2,i+1)+v(4,i+1)]*v(:,i+1);
%     phi=[q(1,i+1)*q(1,i+1)+q(2,i+1)*q(2,i+1)-1;(q(1,i+1)-q(3,i+1))^2+(q(2,i+1)-q(4,i+1))^2-1];%约束
%     phiq=[2*q(1,i+1) 2*q(2,i+1) 0 0;2*(q(1,i+1)-q(3,i+1)) 2*(q(2,i+1)-q(4,i+1)) -2*(q(1,i+1)-q(3,i+1)) -2*(q(2,i+1)-q(4,i+1))];%笛卡尔
%     phiT=phiq*v(:,i+1);
%     P1=P-2*alpha*phiT-(beta^2)*phi;
%     LEFT=[A phiq';phiq zeros(2,2)];
%     RIGHT=[B;P1];
%     X=(LEFT^-1)*RIGHT;
%     a(1,i+1)=X(1);a(2,i+1)=X(2);a(3,i+1)=X(3);a(4,i+1)=X(4);
%     v(:,i+1)=v(:,i)+h*a(:,i+1);
%     q(:,i+1)=q(:,i)+h*v(:,i+1);
end
%% RK-4
for i=1:(length(t)-1)
    %k1
    P=2*[v(1,i) v(2,i) 0 0;v(1,i)-v(3,i) -v(1,i)+v(3,i) v(2,i)-v(4,i) -v(2,i)+v(4,i)]*v(:,i);
    phi=[q(1,i)*q(1,i)+q(2,i)*q(2,i)-1;(q(1,i)-q(3,i))^2+(q(2,i)-q(4,i))^2-1];%约束
    phiq=[2*q(1,i) 2*q(2,i) 0 0;2*(q(1,i)-q(3,i)) 2*(q(2,i)-q(4,i)) -2*(q(1,i)-q(3,i)) -2*(q(2,i)-q(4,i))];%笛卡尔
    phiT=phiq*v(:,i);
    P1=P-2*alpha*phiT-(beta^2)*phi;
    LEFT=[A phiq';phiq zeros(2,2)];
    RIGHT=[B;P1];
    X=(LEFT^-1)*RIGHT;
    ak1(1)=X(1);ak1(2)=X(2);ak1(3)=X(3);ak1(4)=X(4);
    vk1=v(:,i)+h*ak1;
    qk1=q(:,i)+h*vk1;
    %k2
    P=2*[vk1(1) vk1(2) 0 0;vk1(1)-vk1(3) -vk1(1)+vk1(3) vk1(2)-vk1(4) -vk1(2)+vk1(4)]*vk1;
    phi=[qk1(1)*qk1(1)+qk1(2)*qk1(2)-1;(qk1(1)-qk1(3))^2+(qk1(2)-qk1(4))^2-1];%约束
    phiq=[2*qk1(1) 2*qk1(2) 0 0;2*(qk1(1)-qk1(3)) 2*(qk1(2)-qk1(4)) -2*(qk1(1)-qk1(3)) -2*(qk1(2)-qk1(4))];%笛卡尔
    phiT=phiq*vk1;
    P1=P-2*alpha*phiT-(beta^2)*phi;
    LEFT=[A phiq';phiq zeros(2,2)];
    RIGHT=[B;P1];
    X=(LEFT^-1)*RIGHT;
    ak2(1)=X(1);ak2(2)=X(2);ak2(3)=X(3);ak2(4)=X(4);
    vk2=v(:,i)+h*ak2;
    qk2=q(:,i)+h*vk2;
    %k3
    P=2*[vk2(1) vk2(2) 0 0;vk2(1)-vk2(3) -vk2(1)+vk2(3) vk2(2)-vk2(4) -vk2(2)+vk2(4)]*vk2;
    phi=[qk2(1)*qk2(1)+qk2(2)*qk2(2)-1;(qk2(1)-qk2(3))^2+(qk2(2)-qk2(4))^2-1];%约束
    phiq=[2*qk2(1) 2*qk2(2) 0 0;2*(qk2(1)-qk2(3)) 2*(qk2(2)-qk2(4)) -2*(qk2(1)-qk2(3)) -2*(qk2(2)-qk2(4))];%笛卡尔
    phiT=phiq*vk2;
    P1=P-2*alpha*phiT-(beta^2)*phi;
    LEFT=[A phiq';phiq zeros(2,2)];
    RIGHT=[B;P1];
    X=(LEFT^-1)*RIGHT;
    ak3(1)=X(1);ak3(2)=X(2);ak3(3)=X(3);ak3(4)=X(4);
    vk3=v(:,i)+h*ak3;
    qk3=q(:,i)+h*vk3;
    %k4
    P=2*[vk3(1) vk3(2) 0 0;vk3(1)-vk3(3) -vk3(1)+vk3(3) vk3(2)-vk3(4) -vk3(2)+vk3(4)]*vk3;
    phi=[qk3(1)*qk3(1)+qk3(2)*qk3(2)-1;(qk3(1)-qk3(3))^2+(qk3(2)-qk3(4))^2-1];%约束
    phiq=[2*qk3(1) 2*qk3(2) 0 0;2*(qk3(1)-qk3(3)) 2*(qk3(2)-qk3(4)) -2*(qk3(1)-qk3(3)) -2*(qk3(2)-qk3(4))];%笛卡尔
    phiT=phiq*vk3;
    P1=P-2*alpha*phiT-(beta^2)*phi;
    LEFT=[A phiq';phiq zeros(2,2)];
    RIGHT=[B;P1];
    X=(LEFT^-1)*RIGHT;
    ak4(1)=X(1);ak4(2)=X(2);ak4(3)=X(3);ak4(4)=X(4);
    vk4=v(:,i)+h*ak4;
    qk4=q(:,i)+h*vk4;
    %
    a(:,i)=ak1;
    v(:,i+1)=v(:,i)+1/6*h*(ak1+2*ak2+2*ak3+ak4);
    q(:,i+1)=q(:,i)+1/6*h*(vk1+2*vk2+2*vk3+vk4);
    a(:,i+1)=zeros(4,1);
end

%% Figure
subplot(2,2,1),plot(t,q);title('位置');
subplot(2,2,2),plot(t,v);title('速度');
subplot(2,2,3),plot(t,a);title('加速度');
subplot(2,2,4),plot(t,-atan(q(1,:)./q(2,:)));title('角度');