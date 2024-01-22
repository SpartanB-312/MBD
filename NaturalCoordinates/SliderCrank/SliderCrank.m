g = -9.8;
m1 = 1; l1 = 1; J(1) = m1 * (l1 ^ 2)/12;
%m2 = 1; l2 = 2; J2 = m2 * (l2 ^ 2)/3;
m3 = 1;
m2 = 2; l2 = 2; J(2) = (m2-m3) * (l2 ^ 2)/12+m3*((l2/2)^2)/3;
mMat(:,:,1)=diag([m1,m1]);
mMat(:,:,2)=diag([m2,m2]);
q(:,1)=[l1;0;l1+l2];
dq(:,1)=[0;0;0];
D=0.5*[1 0 1 0; 0 1 0 1];
G(1,:)=l1^(-2)*[q(2,1) -q(1,1) -q(2,1) q(1,1)];
G(2,:)=l2^(-2)*[-q(2,1) q(1,1)-q(3,1) q(2,1) -q(1,1)+q(3,1)];
F(:,1)=[0;m1*g];
F(:,2)=[0;m2*g];
for i = 1:2
    A(:,:,i)=D'*mMat(:,:,i)*D+G(i,:)'*J(i)*G(i,:);
    B(:,i)=D'*F(:,i);
end
Phi=[q(1,1)^2-q(2,1)^2-l1^2;
     (q(3,1)-q(1,1))^2+q(2,1)^2-l2^2];
Phiq=2*[q(1,1) q(2,1) 0;
        q(1,1)-q(3,1) q(2,1) q(3,1)-q(1,1)];
dPhi=Phiq*dq;
ALHS=[A(:,:,1),zeros(4);zeros(4),A(:,:,2)];
ALHS=ALHS([3,4,7],[3,4,7]);

h=0.001;T=5;t=0:h:T;
nstep=length(t);
alpha = 25; beta = 50; %经验值5~50



for i=1:nstep-1

    D=0.5*[1 0 1 0; 0 1 0 1];
    G(1,:)=l1^(-2)*[q(2,i) -q(1,i) -q(2,i) q(1,i)];
    G(2,:)=l2^(-2)*[-q(2,i) q(1,i)-q(3,i) q(2,i) -q(1,i)+q(3,i)];

    for j = 1:2
        A(:,:,j)=D'*mMat(:,:,j)*D+G(j,:)'*J(j)*G(j,:);
        B(:,j)=D'*F(:,j);
    end
    Phi=[q(1,i)^2+q(2,i)^2-l1^2;
         (q(3,i)-q(1,i))^2+q(2,i)^2-l2^2];
    Phiq=2*[q(1,i) q(2,i) 0;
            q(1,i)-q(3,i) q(2,i) q(3,i)-q(1,i)];
    dPhi=Phiq*dq(:,i);
    P=-2*[dq(1,i),dq(2,i),0;dq(1,i)-dq(3,i),dq(2,i),dq(3,i)-dq(1,i)]*dq(:,i);
    P1 = P - 2 * alpha * dPhi - beta ^ 2 * Phi; %稳定形式
    ALHS=[A(:,:,1),zeros(4);zeros(4),A(:,:,2)];
    temp1=ALHS([7],[5,6]);
    temp2=ALHS([5,6],[7]);
%     ALHS=ALHS([3,4,7],[3,4,7]);
    ALHS=ALHS([5,6,7],[5,6,7]);
%     ALHS([3],[1,2])=temp1;
%     ALHS([1,2],[3])=temp2;
    BRHS=[B(:,1);B(:,2)];
    BRHS=BRHS([5,6,7]);
    LHS=[ALHS,Phiq';Phiq,zeros(2)];
    RHS=[BRHS;P1];
    temp=LHS\RHS;
    ddq(1:3,i)=temp(1:3);
    dq(:,i+1)=dq(:,i)+ddq(:,i)*h;
    q(:,i+1)=q(:,i)+dq(:,i)*h;
end
%% figure
plot(t,q(3,:))