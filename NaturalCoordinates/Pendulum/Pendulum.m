l=1;
m1=1;J1=[1/3*m1*(2*l)^2];
M1=[m1,0;0,m1];
g=-9.8;
q(:,1)=[0.5; -sqrt(3) / 2;0;0];
dq(:,1)=[0;0;0;0];
D=[1 0 0 0;
   0 1 0 0];
G=l^(-2) * [-q(2,1) q(1,1) q(2,1) -q(1,1)];
Phi=q(1,1)^2+q(2,1)^2-l^2;
Phiq=2*[q(1,1) q(2,1) 0 0];
P=-2*(dq(1,1)^2+dq(2,1)^2);

h=0.001;T=5;t=0:h:T;
nstep=length(t);
alpha = 25; beta = 50; %经验值5~50
for i=1:nstep-1
    F=[0;m1*g];
    M=[m1*g*q(1,i)];

    D=[1 0 0 0;
       0 1 0 0];
    G=l^(-2) * [q(4,i)-q(2,i) q(1,i)-q(3,i) q(2,i)-q(4,i) q(3,i)-q(1,i)];
    dG=l^(-2) * [dq(4,i)-dq(2,i) dq(1,i)-dq(3,i) dq(2,i)-dq(4,i) dq(3,i)-dq(1,i)];
%     Phi=[q(3,i)^2+q(4,i)^2;
%          (q(1,i)-q(3,i))^2+(q(2,i)-q(4,i))^2-l^2];
%     Phiq=2*[0 0 q(3,i) q(4,i);
%             (q(1,i)-q(3,i)) (q(2,i)-q(4,i)) -(q(1,i)-q(3,i)) -(q(2,i)-q(4,i))];
%     Phit=[0;0];
    Phi=q(1,i)^2+q(2,i)^2-l^2;
    Phiq=[q(1,i) q(2,i) 0 0];
    
    Phit=0;
    PhiT=Phiq*dq(:,i)+Phit;
%     P=[0;-2*(dq(1,i)^2+dq(2,i)^2)];
    P=-2*(dq(1,i)^2+dq(2,i)^2);
    P1 = P - 2 * alpha * PhiT - beta ^ 2 * Phi; %稳定形式
    A=D'*M1*D+G'*J1*G;
    B=D'*F+G'*M;
    Phiq=[q(1,i) q(2,i)];
    LHS=[A(1:2,1:2),Phiq';Phiq,zeros(1)];
    RHS=[B(1:2);P1];
    temp=LHS\RHS;
    ddq(1:2,i)=temp(1:2);
    ddq(3:4,i)=zeros(2,1);
    dq(:,i+1)=dq(:,i)+ddq(:,i)*h;
    q(:,i+1)=q(:,i)+dq(:,i)*h;
end
%% figure
plot(q(1,:))