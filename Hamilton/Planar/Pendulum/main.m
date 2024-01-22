%
clear;
m1=1;J1=0;
M = [m1 0 0;
    0 m1 0;
    0 0 J1];
g=-9.81;
q(:,1)=[0.5; -sqrt(3) / 2; pi / 6];
dq(:,1)=[0;0;0];
sigma(:,1)=[0;0];
Phiq=PhiqCal(q(:,1));
pstar(:,1)=M*dq(:,1)+Phiq'*sigma(:,1);
%
Time=5;h=0.001;nstep=Time/h;
alpha = 25; beta = 25; %经验值5~50
%
tic
for i=1:nstep-1
    %
    Phiq=PhiqCal(q(:,i));
    Phi(:,i)=PhiCal(q(:,i));
    Qe=[0;m1*g;0];
    LHS=[M,Phiq';
        Phiq,zeros(2)];
    RHS=[pstar(:,i);zeros(2,1)];
    temp=LHS\RHS;
    dq(:,i)=temp(1:3);
    sigma(:,i)=temp(4:5);
    %
    dPhiq=dPhiqCal(q(:,i),dq(:,i));
    %
    dpstar(:,i)=Qe+dPhiq'*sigma(:,i);
    %
    q(:,i+1)=q(:,i)+h*dq(:,i);
    pstar(:,i+1)=pstar(:,i)+h*dpstar(:,i);
    %
end
Phi(:,i+1)=PhiCal(q(:,i+1));
toc
%
%% Figure
t=h:h:Time;
%
figure(1)
hold on
grid on
plot(t,q(3,:),'-g','LineWidth',1);
hold off;
%
figure(2)
hold on
grid on
plot(t,Phi(1,:),'-g','LineWidth',1);
plot(t,Phi(2,:),'-r','LineWidth',1);
legend('Phi1','Phi2')
hold off