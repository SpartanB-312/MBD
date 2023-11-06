l=1;A=0.02*0.02;
E=6*10^7;nu=0.3;
rho=6000;
m=rho*A*l;g=-9.81;

mu=E/(2*(1+nu));
lam=nu*E/((1-2*nu)*(1+nu));

Time = 0.01; h = 0.001; nstep = Time / h;
alpha=25;beta=50;

e=[0,0,0,-1,0,-1,0,-1]';
%e=[0,0,1,0,1,0,1,0]';
de=[0,0,0,0,0,0,0,0]';
M=MassCal(rho,l);
for i = 1:nstep - 1
    K=KCal(e(:,i),lam,mu,A,l);
    Fe=-K*e(:,i);
    Qg = QgCal(rho,g,A,l);

    Phi=[e(1,i);
         e(2,i)];
    Phie=[1,0,0,0,0,0,0,0;
          0,1,0,0,0,0,0,0];
    dPhi=Phie*de(:,i);
    dPhie=[0,0,0,0,0,0,0,0;
           0,0,0,0,0,0,0,0];

    P=-dPhie*de(:,i);
    P1=P-2*alpha*Phi-2*beta*dPhi;

    LHS=[M,Phie';Phie,zeros(2)];
    RHS=[Qg-Fe;P1];
    temp=LHS\RHS;
    dde(:,i)=temp(1:8);
    de(:,i+1)=de(:,i)+h*dde(:,i);
    e(:,i+1)=e(:,i)+h*de(:,i);
end