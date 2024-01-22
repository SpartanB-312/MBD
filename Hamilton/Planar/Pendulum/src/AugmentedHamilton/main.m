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
Time=5;h=0.01;nstep=Time/h;
alpha = 1e9; omega = 0.1; ksi = 1000;
%
tic
for i=1:nstep-1
    %
    if i > 1
        dq(:,i)=dq(:,i-1);
        sigma(:,i)=sigma(:,i-1);
    end
    %dq(:,i)=zeros(3,1);
    %sigma(:,i)=zeros(2,1);
    %q(:,i)=q(:,i-1);
    %pstar(:,i)=pstar(:,i-1);
    j=1;
    while max(abs(1)) > 1e-5
        Phiq=PhiqCal(q(:,i));
        dPhi=Phiq*dq(:,i);
        dPhiq=dPhiqCal(q(:,i),dq(:,i));
        Phi(:,i)=PhiCal(q(:,i));
        sigma(:,i)=sigma(:,i)+alpha*(dPhi+2*ksi*omega*Phi(:,i)+(omega^2)*h*sum(Phi,2));
        
        Qe=[0;m1*g;0];
        LHS=M+Phiq'*alpha*Phiq;
        penalty=Phiq'*alpha*(2*ksi*omega*Phi(:,i)+(omega^2)*h*sum(Phi,2))+Phiq'*sigma(:,i);
        RHS=pstar(:,i)-penalty;
        dq(:,i)=LHS\RHS;
        dpstar(:,i)=Qe+dPhiq'*alpha*(dPhi+2*ksi*omega*Phi(:,i)+(omega^2)*h*sum(Phi,2))+dPhiq'*sigma(:,i);
        
        q(:,i+1)=q(:,i)+h*dq(:,i);
        pstar(:,i+1)=pstar(:,i)+h*dpstar(:,i);
        
        if j > 5
            break;
        end
        j = j+1;
    end
    
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