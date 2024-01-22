%4.9
m1=1;M(:,:,1)=diag([1,1,1/4]);
m2=1;M(:,:,2)=diag([1,1,1/4]);
g=-9.8;
m1=1;M(:,:,1)=diag([1,1,0]);
m2=1;M(:,:,2)=diag([1,1,0]);

I2=eye(2);I3=eye(3);zero3=zeros(3);
I2tilde=[0 -1;1 0];

T=[I3,zero3;I3,I3];


ai(:,1)=[0;0];
ai(:,2)=[0;0];
bi(:,1)=[1;0];
bi(:,2)=[1;0];

q(:,1,1)=[0.5;-sqrt(3)/2;-pi/3];
q(:,2,1)=[1;-sqrt(3);-pi/3];
z(:,1)=[-pi/3;0];
A(:,:,1)=CgaCal(q(3,1,1));
A(:,:,2)=CgaCal(q(3,2,1));
Qe(1:3,1)=[0;m1*g;0];
Qe(4:6,1)=[0;m2*g;0];

dz(:,1)=[0;0];

rj(:,1)=[0;0];
rj(:,2)=q(1:2,1,1)+A(:,:,1)*ai(:,1);



B1_0=[I2 I2tilde*(q(1:2,1,1));zeros(1,2) 1];
B2_1=[I2 I2tilde*(q(1:2,2,1)-q(1:2,1,1));zeros(1,2) 1];
V1_0=[I2tilde*(q(1:2,1,1)-rj(:,1));1];
V2_0=[I2tilde*(q(1:2,2,1)-rj(:,1));1];
V2_1=[I2tilde*(q(1:2,2,1)-rj(:,2));1];
V=[V1_0 zeros(3,1);V2_0 V2_1];
temp=V*dz;
dq(:,1,1)=temp(1:3);
dq(:,2,1)=temp(4:6);

Time=5;h=0.001;nstep=Time/h;
for i=1:nstep-1

    A(:,:,1)=CgaCal(q(3,1,i));
    A(:,:,2)=CgaCal(q(3,2,i));
    rj(:,1)=[0;0];
    rj(:,2)=q(1:2,1,i)+A(:,:,1)*ai(:,1);
    
    V1_0=[I2tilde*(q(1:2,1,i)-rj(:,1));1];
    V2_0=[I2tilde*(q(1:2,2,i)-rj(:,1));1];
    V2_1=[I2tilde*(q(1:2,2,i)-rj(:,2));1];
    B1_0=[I2 I2tilde*(q(1:2,1,i));zeros(1,2) 1];
    B2_1=[I2 I2tilde*(q(1:2,2,i)-q(1:2,1,i));zeros(1,2) 1];
    V=[V1_0 zeros(3,1);V2_0 V2_1];
    
    temp=V*dz(:,i);
    dq(:,1,i)=temp(1:3);
    dq(:,2,i)=temp(4:6);
    drj(:,1)=[0;0];
    drj(:,2)=dq(1:2,1,i)+(I2tilde.*dq(3,1,i))*A(:,:,1)*ai(:,1);
    dV1_0=[I2tilde*(dq(1:2,1,i)-drj(:,1));0];
    dV2_0=[I2tilde*(dq(1:2,2,i)-drj(:,1));0];
    dV2_1=[I2tilde*(dq(1:2,2,i)-drj(:,2));0];
    
    dV=[dV1_0 zeros(3,1);dV2_0 dV2_1];
    Mbar=V'*[M(:,:,1),zeros(3);zeros(3),M(:,:,2)]*V;
    Qebar=V'*Qe;
    Qvbar=-V'*[M(:,:,1),zeros(3);zeros(3),M(:,:,2)]*dV*dz(:,i);

    ddz(:,i)=(Mbar^-1)*(Qebar+Qvbar);
    dz(:,i+1)=dz(:,i)+h*ddz(:,i);
    z(:,i+1)=z(:,i)+h*dz(:,i);

    q(:,1,i+1)=[CgaCal(z(1,i+1))*bi(:,1);z(1,i+1)];
    q(:,2,i+1)=q(:,1,i+1)+[CgaCal(q(3,1,i+1))*ai(:,1)+CgaCal(q(3,1,i+1)+z(2,i+1))*bi(:,2);z(2,i+1)];
end
%% figure
figure(1)
theta14fig=reshape(q(3,1,:),1,nstep);
theta24fig=reshape(q(3,2,:),1,nstep);
hold on
grid on
plot(90+theta14fig*180/pi,'-g');
plot(90+theta24fig*180/pi,'-r');
legend('theta1','theta2')
hold off
figure(2)
load('DoublePendulum.mat');
plot((pi/2+theta14fig)*180/pi-MEA_ANGLE_1Q);