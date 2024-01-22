clear all;
l=1;h=0.02;w=0.02;
E=6*10^7;nu=0.3;
rho=6000;
m=rho*h*l*w;g=-9.81;

mu=E/(2*(1+nu));
lam=nu*E/((1-2*nu)*(1+nu));

Time = 0.05; Timeh = 0.001; nstep = Time / Timeh;
alpha=25;beta=50;

e=[0,0,1,0,0,-l/2,1,0,0,-l,1,0]';
%e=[0,0,0,1,l/2,0,0,1,l,0,0,1]';
%e=[0,0,1,0,1,0,1,0]';
de=[0,0,0,0,0,0,0,0,0,0,0,0]';
M=MassCal(rho,h,l,w);
for i = 1:nstep - 1
    K=KCal(e(:,i),lam,mu,h,l,w);
    Fe=-K*e(:,i);
    Qg = QgCal(rho,g,h,l,w);

    Phi=[e(1,i);
         e(2,i)];
    Phie=[1,0,0,0,0,0,0,0,0,0,0,0;
          0,1,0,0,0,0,0,0,0,0,0,0];
    dPhi=Phie*de(:,i);
    dPhie=[0,0,0,0,0,0,0,0,0,0,0,0;
           0,0,0,0,0,0,0,0,0,0,0,0];

    P=-dPhie*de(:,i);
    P1=P-2*alpha*Phi-2*beta*dPhi;

    LHS=[M,Phie';Phie,zeros(2)];
    RHS=[Qg-Fe;P1];
    temp=LHS\RHS;
    dde(:,i)=temp(1:12);
    de(:,i+1)=de(:,i)+Timeh*dde(:,i);
    e(:,i+1)=e(:,i)+Timeh*de(:,i);
end