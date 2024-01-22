%
%约束 Phiz=Phiq*(dq/dz)
m1=1;M(:,:,1)=diag([1,1,1/12*m1]);
m2=1;M(:,:,2)=diag([1,1,1/12*m2*4]);
m3=1;M(:,:,3)=diag([1,1,1]);
g=-9.8;

I2=eye(2);I3=eye(3);zero3=zeros(3);
I2tilde=[0 -1;1 0];

T=[I3,zero3,zero3;I3,I3,zero3;I3,I3,I3];

ai(:,1)=[0.5;0];
ai(:,2)=[1;0];
ai(:,3)=[0;0];
bi(:,1)=[0.5;0];
bi(:,2)=[1;0];
bi(:,3)=[0;0];


q(:,1,1)=[0.5;0;0];
q(:,2,1)=[2;0;0];
q(:,3,1)=[3;0;0];
z(:,1)=[0;0;0];
A(:,:,1)=CgaCal(q(3,1,1));
A(:,:,2)=CgaCal(q(3,2,1));
A(:,:,3)=CgaCal(q(3,3,1));
Qe(1:3,1)=[0;m1*g;0];
Qe(4:6,1)=[0;m2*g;0];
Qe(7:9,1)=[0;m3*g;0];

dz(:,1)=[0;0;0];

rj(:,1)=[0;0];
rj(:,2)=q(1:2,1,1)+A(:,:,1)*ai(:,1);
rj(:,3)=q(1:2,2,1)+A(:,:,2)*ai(:,2);

%Bi_i-1
B1_0=[I2 I2tilde*(q(1:2,1,1));zeros(1,2) 1];
B2_1=[I2 I2tilde*(q(1:2,2,1)-q(1:2,1,1));zeros(1,2) 1];
B3_2=[I2 I2tilde*(q(1:2,3,1)-q(1:2,2,1));zeros(1,2) 1];
%Vi_j
V1_0=[I2tilde*(q(1:2,1,1)-rj(:,1));1];
V2_0=[I2tilde*(q(1:2,2,1)-rj(:,1));1];
V2_1=[I2tilde*(q(1:2,2,1)-rj(:,2));1];
V3_0=[I2tilde*(q(1:2,3,1)-rj(:,1));1];
V3_1=[I2tilde*(q(1:2,3,1)-rj(:,2));1];
V3_2=[I2tilde*(q(1:2,3,1)-rj(:,3));1];
V=[V1_0,zeros(3,1),zeros(3,1);V2_0,V2_1,zeros(3,1);V3_0,V3_1,V3_2];

temp=V*dz;
dq(:,1,1)=temp(1:3);
dq(:,2,1)=temp(4:6);
dq(:,3,1)=temp(7:9);

Phi=[q(2,3,1);q(3,3,1)];
Phiq=[0,0,0,0,0,0,0,1,0;
      0,0,0,0,0,0,0,0,1];
Phiz=Phiq*V;
dPhiq=zeros(6,2);

Time=5;h=0.001;nstep=Time/h;
alpha=25;beta=50;
for i=1:nstep-1

    A(:,:,1)=CgaCal(q(3,1,i));
    A(:,:,2)=CgaCal(q(3,2,i));
    A(:,:,3)=CgaCal(q(3,3,i));
    rj(:,1)=[0;0];
    rj(:,2)=q(1:2,1,i)+A(:,:,1)*ai(:,1);
    rj(:,3)=q(1:2,2,i)+A(:,:,2)*ai(:,2);
    
    V1_0=[I2tilde*(q(1:2,1,i)-rj(:,1));1];
    V2_0=[I2tilde*(q(1:2,2,i)-rj(:,1));1];
    V2_1=[I2tilde*(q(1:2,2,i)-rj(:,2));1];
    V3_0=[I2tilde*(q(1:2,3,i)-rj(:,1));1];
    V3_1=[I2tilde*(q(1:2,3,i)-rj(:,2));1];
    V3_2=[I2tilde*(q(1:2,3,i)-rj(:,3));1];
    B1_0=[I2 I2tilde*(q(1:2,1,i));zeros(1,2) 1];
    B2_1=[I2 I2tilde*(q(1:2,2,i)-q(1:2,1,i));zeros(1,2) 1];
    B3_2=[I2 I2tilde*(q(1:2,3,i)-q(1:2,2,i));zeros(1,2) 1];
    V=[V1_0,zeros(3,1),zeros(3,1);V2_0,V2_1,zeros(3,1);V3_0,V3_1,V3_2];
    
    temp=V*dz(:,i);
    dq(:,1,i)=temp(1:3);
    dq(:,2,i)=temp(4:6);
    dq(:,3,1)=temp(7:9);
    drj(:,1)=[0;0];
    drj(:,2)=dq(1:2,1,i)+(I2tilde.*dq(3,1,i))*A(:,:,1)*ai(:,1);
    drj(:,3)=dq(1:2,2,i)+(I2tilde.*dq(3,2,i))*A(:,:,2)*ai(:,2);
    dV1_0=[I2tilde*(dq(1:2,1,i)-drj(:,1));0];
    dV2_0=[I2tilde*(dq(1:2,2,i)-drj(:,1));0];
    dV2_1=[I2tilde*(dq(1:2,2,i)-drj(:,2));0];
    dV3_0=[I2tilde*(dq(1:2,3,i)-drj(:,1));0];
    dV3_1=[I2tilde*(dq(1:2,3,i)-drj(:,2));0];
    dV3_2=[I2tilde*(dq(1:2,3,i)-drj(:,3));0];
    
    dV=[dV1_0,zeros(3,1),zeros(3,1);dV2_0,dV2_1,zeros(3,1);dV3_0,dV3_1,dV3_2];

    Phi=[q(2,3,i);q(3,3,i)];
    Phiq=[0,0,0,0,0,0,0,1,0;
          0,0,0,0,0,0,0,0,1];
    Phiz=Phiq*V;
    dPhiq=zeros(2,9);
    dPhiz=dPhiq*V+Phiq*dV;
    Phit=[0;0];
    dPhit=[0;0];
    dPhi=Phiz*dz(:,i);
    
    P=-(dPhiz*dz(:,i)+dPhit);
    P1=P-2*alpha*Phi-2*beta*dPhi;
    
    Mbar=V'*[M(:,:,1),zero3,zero3;zero3,M(:,:,2),zero3;zero3,zero3,M(:,:,3)]*V;
    Qebar=V'*Qe;
    Qvbar=-V'*[M(:,:,1),zero3,zero3;zero3,M(:,:,2),zero3;zero3,zero3,M(:,:,3)]*dV*dz(:,i);

    LHS=[Mbar,Phiz';Phiz,zeros(2)];
    RHS=[Qebar+Qvbar;P1];

    %ddz(:,i)=(Mbar^-1)*(Qebar+Qvbar);
    temp=LHS\RHS;
    ddz(:,i)=temp(1:3);
    dz(:,i+1)=dz(:,i)+h*ddz(:,i);
    z(:,i+1)=z(:,i)+h*dz(:,i);

    q(:,1,i+1)=[CgaCal(z(1,i+1))*bi(:,1);z(1,i+1)];
    q(:,2,i+1)=q(:,1,i+1)+[CgaCal(q(3,1,i+1))*ai(:,1)+CgaCal(q(3,1,i+1)+z(2,i+1))*bi(:,2);z(2,i+1)];
    q(:,3,i+1)=q(:,2,i+1)+[CgaCal(q(3,2,i+1))*ai(:,2)+CgaCal(q(3,2,i+1)+z(3,i+1))*bi(:,3);z(3,i+1)];
end
%% figure
figure(1)
x34fig=reshape(q(1,3,:),1,nstep);
hold on
grid on
plot(x34fig,'-g');
hold off