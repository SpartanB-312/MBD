%%
h = 0.001; T = 5; t = 0:h:T-h; %步长及时间
deltal=0.1;l0=1;nl=11;
deltam=0.1;m0=1;nm=11;
data4save=zeros(5000,nl*nm);
for i = 1:nl
    for j = 1:nm
        angle2t = Pendulum(m0+deltam*(j-1),l0+deltal*(i-1));
        data4save(:,nl*(i-1)+j) = angle2t'*180/pi;
    end
end
%%
%writematrix(data4save,'PendulumData.dat','Delimiter',',');
csvwrite('PendulumData.csv',data4save);