%4.9
m1 = 1; M(:, :, 1) = diag([1, 1, 1/4]);
m2 = 1; M(:, :, 2) = diag([1, 1, 1/4]);
m3 = 1; M(:, :, 3) = diag([1, 1, 1]);
g = -9.8;


I2 = eye(2); I3 = eye(3); zero3 = zeros(3);zero2=zeros(2);
I2tilde = [0 -1; 1 0];

T=[I3,zero3,zero3;I3,I3,zero3;I3,I3,I3];

r1_1 = [-0.5; 0];
r1_2 = [0.5; 0];
r2_2 = [-1; 0];
r2_3 = [1; 0];
r3_3 = [0;0];

q(:, 1, 1) = [0.5;0;0];
q(:, 2, 1) = [2;0;0];
q(:, 3, 1) = [3;0;0];
A(:, :, 1) = CgaCal(q(3, 1, 1));
A(:, :, 2) = CgaCal(q(3, 2, 1));
A(:, :, 3) = CgaCal(q(3, 3, 1));
rj(:, 1) = [0; 0];
rj(:, 2) = q(1:2, 1, 1) + A(:, :, 1) * r1_2(:, 1);
rj(:,3) = q(1:2, 2, 1) + A(:, :, 2) * r2_3(:, 1);
bj(:, 1) = [-I2tilde * rj(:, 1); 1];
bj(:, 2) = [-I2tilde * rj(:, 2); 1];
bj(:, 3) = [-I2tilde * rj(:, 3); 1];
Rd = [bj(:, 1), zeros(3, 1), zeros(3,1); zeros(3, 1), bj(:, 2), zeros(3,1); zeros(3,1), zeros(3,1), bj(:,3)];

%Z(:, 1, 1) = [CgaCal(z(1,1)) * r1_1;z(1,1)];
%Z(:, 2, 1) = Z(:, 1, 1) + [(CgaCal(z(1,1)) * r1_2 - CgaCal(Z(3,1,1)+z(2,1)) * r2_2);z(2,1)];

dz(:, 1) = [0; 0; 0];
dj(:, 1, 1) = [(dz(1, 1) ^ 2) * rj(:, 1); 0];
dj(:, 2, 1) = [(dz(2, 1) ^ 2 - dz(1, 1) ^ 2) * rj(:, 2); 0];
dj(:, 3, 1) = [(dz(3, 1) ^ 2 - dz(2, 1) ^ 2) * rj(:, 2); 0];
temp = T * Rd * dz(:, 1);
dZ(:, 1, 1) = temp(1:3);
dZ(:, 2, 1) = temp(4:6);
dZ(:, 3, 1) = temp(7:9);
%temp = T * (Rd * ddz(:, 1) + [dj(:, 1, 1); dj(:, 2, 1)]);
%ddZ(:, 1, 1) = temp(1:3);
%ddZ(:, 2, 1) = temp(4:6);


ei(:, 1, 1) = [-dZ(3, 1, 1) ^ 2 * q(1:2, 1, 1); 0];
ei(:, 2, 1) = [-dZ(3, 2, 1) ^ 2 * q(1:2, 2, 1); 0];
ei(:, 2, 1) = [-dZ(3, 3, 1) ^ 2 * q(1:2, 3, 1); 0];

Di(:, :, 1) = [I2, I2tilde * q(1:2, 1, 1); zeros(1, 2), 1];
Di(:, :, 2) = [I2, I2tilde * q(1:2, 2, 1); zeros(1, 2), 1];
Di(:, :, 3) = [I2, I2tilde * q(1:2, 3, 1); zeros(1, 2), 1];
D = [Di(:, :, 1), zero3, zero3; zero3, Di(:, :, 2), zero3; zero3, zero3, Di(:, :, 3)];
dq(:, 1, 1) = Di(:, :, 1) * dZ(:, 1, 1);
dq(:, 2, 1) = Di(:, :, 2) * dZ(:, 2, 1);
dq(:, 3, 1) = Di(:, :, 3) * dZ(:, 3, 1);

z(:, 1) = [0; 0; 0];

Qe(1:3, 1) = [0; m1 * g; 0];
Qe(4:6, 1) = [0; m2 * g; 0];
Qe(7:9, 1) = [0; m3 * g; 0];

Time = 5; h = 0.001; nstep = Time / h;
alpha=25;beta=50;
for i = 1:nstep - 1

    A(:, :, 1) = CgaCal(q(3, 1, i));
    A(:, :, 2) = CgaCal(q(3, 2, i));
    A(:, :, 3) = CgaCal(q(3, 3, i));
    rj(:, 1) = [0; 0];
    rj(:, 2) = q(1:2, 1, i) + A(:, :, 1) * r1_2;
    rj(:,3) = q(1:2, 2, i) + A(:, :, 2) * r2_3;
    bj(:, 1) = [-I2tilde * rj(:, 1); 1];
    bj(:, 2) = [-I2tilde * rj(:, 2); 1];
    bj(:, 3) = [-I2tilde * rj(:, 3); 1];
    Rd = [bj(:, 1), zeros(3, 1), zeros(3,1); zeros(3, 1), bj(:, 2), zeros(3,1); zeros(3,1), zeros(3,1), bj(:,3)];

    temp = T * Rd * dz(:, i);
    dZ(:, 1, i) = temp(1:3);
    dZ(:, 2, i) = temp(4:6);
    dZ(:, 3, i) = temp(7:9);
    
    dj(:, 1, i) = [(dZ(3,1, i) ^ 2) * rj(:, 1); 0];
    dj(:, 2, i) = [(dZ(3,2, i) ^ 2 - dZ(3,1, i) ^ 2) * rj(:, 2); 0];
    dj(:, 3, i) = [(dZ(3,3, i) ^ 2 - dZ(3,2, i) ^ 2) * rj(:, 3); 0];

    ei(:, 1, i) = [-dZ(3, 1, i) ^ 2 * q(1:2, 1, i); 0];
    ei(:, 2, i) = [-dZ(3, 2, i) ^ 2 * q(1:2, 2, i); 0];
    ei(:, 3, i) = [-dZ(3, 3, i) ^ 2 * q(1:2, 3, i); 0];

    Di(:, :, 1) = [I2, I2tilde * q(1:2, 1, i); zeros(1, 2), 1];
    Di(:, :, 2) = [I2, I2tilde * q(1:2, 2, i); zeros(1, 2), 1];
    Di(:, :, 3) = [I2, I2tilde * q(1:2, 3, i); zeros(1, 2), 1];
    D = [Di(:, :, 1), zero3, zero3; zero3, Di(:, :, 2), zero3; zero3, zero3, Di(:, :, 3)];

    dq(:, 1, i) = Di(:, :, 1) * dZ(:, 1, i);
    dq(:, 2, i) = Di(:, :, 2) * dZ(:, 2, i);
    dq(:, 3, i) = Di(:, :, 3) * dZ(:, 3, i);

    Phi=[q(2,3,i);q(3,3,i)];
    Phiq=[0,0,0,0,0,0,0,1,0;
          0,0,0,0,0,0,0,0,1];
    Phiz=Phiq*(D*T*Rd);
    dPhiq=zeros(2,9);
    dDi(:, :, 1) = [zero2, I2tilde * dq(1:2, 1, i); zeros(1, 2), 0];
    dDi(:, :, 2) = [zero2, I2tilde * dq(1:2, 2, i); zeros(1, 2), 0];
    dDi(:, :, 3) = [zero2, I2tilde * dq(1:2, 3, i); zeros(1, 2), 0];
    dD = [dDi(:, :, 1), zero3, zero3; zero3, dDi(:, :, 2), zero3; zero3, zero3, dDi(:, :, 3)];
    drj(:, 1) = [0; 0];
    drj(:, 2) = dq(1:2, 1, i) + A(:, :, 1) * r1_2;
    drj(:, 3) = dq(1:2, 2, i) + A(:, :, 2) * r2_3;
    dbj(:, 1) = [-I2tilde * drj(:, 1); 1];
    dbj(:, 2) = [-I2tilde * drj(:, 2); 1];
    dbj(:, 3) = [-I2tilde * drj(:, 3); 1];
    dRd = [dbj(:, 1), zeros(3, 1), zeros(3,1);
           zeros(3, 1), dbj(:, 2), zeros(3,1);
           zeros(3,1), zeros(3,1), dbj(:, 3)];
    dDTRd=dD*T*Rd+D*T*dRd;
    dPhiz=dPhiq*(D*T*Rd)+Phiq*dDTRd;
    Phit=[0;0];
    dPhit=[0;0];
    dPhi=Phiz*dz(:,i);
    
    P=-(dPhiz*dz(:,i)+dPhit);
    P1=P-2*alpha*Phi-2*beta*dPhi;

    Mbar = Rd' * T' * D' * [M(:, :, 1), zero3, zero3; zero3, M(:, :, 2), zero3;zero3,zero3,M(:, :, 3)] * D * T * Rd;
    Qebar = Rd' * T' * D' * Qe;
    Qvbar = Rd' * T' * D' * (- [M(:, :, 1), zero3, zero3; zero3, M(:, :, 2), zero3;zero3,zero3,M(:, :, 3)] * ([ei(:, 1, i); ei(:, 2, i); ei(:, 3, i)] + D * T * [dj(:, 1, i); dj(:, 2, i); dj(:, 3, i)]));

    LHS=[Mbar,Phiz';Phiz,zeros(2)];
    RHS=[Qebar+Qvbar;P1];
    %ddz(:, i) = (Mbar ^ -1) * (Qebar + Qvbar);
    temp=LHS\RHS;
    ddz(:,i)=temp(1:3);
    dz(:, i + 1) = dz(:, i) + h * ddz(:, i);
    z(:, i + 1) = z(:, i) + h * dz(:, i);

    q(:, 1, i+1) = [-CgaCal(z(1,i+1)) * r1_1;z(1,i+1)];
    q(:, 2, i+1) = q(:, 1, i+1) + [(CgaCal(z(1,i+1)) * r1_2 - CgaCal(q(3,1,i+1)+z(2,i+1)) * r2_2);z(2,i+1)];
    q(:, 3, i+1) = q(:, 2, i+1) + [(CgaCal(q(3,2,i+1)) * r2_3 - CgaCal(q(3,2,i+1)+z(3,i+1)) * r3_3);z(3,i+1)];

end

%% figure
figure(1)
x34fig=reshape(q(1,3,:),1,nstep);
hold on
grid on
plot(x34fig,'-g');
hold off
