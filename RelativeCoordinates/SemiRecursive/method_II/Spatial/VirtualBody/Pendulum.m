%4.9
clear;
g = -9.8;
m1 = 0; M(:, :, 1) = diag([m1, m1, m1, 0, 0, 0]);
m2 = 1; M(:, :, 2) = diag([m2, m2, m2, 0, 0, 0]);

I2 = eye(2); I3 = eye(3); I6 = eye(6);
zero3 = zeros(3); zero6 = zeros(6);
I2tilde = [0 -1; 1 0];

T = [I6, zero6; I6, I6];

v_L=0;%虚拟长度
r1_1 = [-v_L; 0; 0];
r1_2 = [0; 0; 0];
r2_2 = [-1; 0; 0];
nj1_1 = [0; 1; 0];
nj1_2 = [0; 0; 1];
nj2_2 = [0; 0; 1];

z(:, 1) = [-pi / 3; 0];
q(:, 1, 1) = [0.5*v_L; -sqrt(3) / 2*v_L; 0; 0; 0; -pi / 3];
q(:, 2, 1) = [1/2+v_L*0.5; -sqrt(3)/2*(1+v_L); 0; 0; 0;-pi / 3];

A(:, :, 1) = CgaCal(q(4, 1, 1),q(5, 1, 1),q(6, 1, 1));
A(:, :, 2) = CgaCal(q(4, 2, 1),q(5, 2, 1),q(6, 2, 1));
nj_1=A(:,:,1)*nj1_1;
nj_2=A(:,:,2)*nj2_2;
rj(:, 1) = [0; 0; 0];
rj(:, 2) = q(1:3, 1, 1) + A(:, :, 1) * r1_2(:, 1);
bj(:, 1) = [Vec2Mat(rj(:, 1))*nj_1; nj_1];
bj(:, 2) = [Vec2Mat(rj(:, 2))*nj_2; nj_2];
Rd = [bj(:, 1), zeros(6, 1); zeros(6, 1), bj(:, 2)];


dz(:, 1) = [0; 0];
temp = T * Rd * dz(:, 1);
dZ(:, 1, 1) = temp(1:6);
dZ(:, 2, 1) = temp(7:12);
wi(:,1)=dZ(4:6, 1, 1);
wi(:,2)=dZ(4:6, 2, 1);
wimat(:,:,1)=Vec2Mat(wi(:,1));
wimat(:,:,2)=Vec2Mat(wi(:,2));

dj(:, 1, 1) = [(-wimat(:,:,1)*wimat(:,:,1)) * rj(:, 1); wimat(:,:,1)*nj_1*dz(1,1)];
dj(:, 2, 1) = [(wimat(:,:,1)*wimat(:,:,1)-wimat(:,:,2)*wimat(:,:,2)) * rj(:, 2); wimat(:,:,2)*nj_2*dz(2,1)];


ei(:, 1, 1) = [wimat(:,:,1) ^ 2 * q(1:3, 1, 1); zeros(3,1)];
ei(:, 2, 1) = [wimat(:,:,2) ^ 2 * q(1:3, 2, 1); zeros(3,1)];

Di(:, :, 1) = [I3, -Vec2Mat(q(1:3, 1, 1)); zeros(3), I3];
Di(:, :, 2) = [I3, -Vec2Mat(q(1:3, 2, 1)); zeros(3), I3];
D = [Di(:, :, 1), zero6; zero6, Di(:, :, 2)];
dq(:, 1, 1) = Di(:, :, 1) * dZ(:, 1, 1);
dq(:, 2, 1) = Di(:, :, 1) * dZ(:, 2, 1);

Qe(1:6, 1) = [0; m1 * g; 0; 0; 0; 0];
Qe(7:12, 1) = [0; m2 * g; 0; 0; 0; 0];

Time = 5; h = 0.001; nstep = Time / h;

for i = 1:nstep - 1

    A(:, :, 1) = CgaCal(q(4, 1, i),q(5, 1, i),q(6, 1, i));
    A(:, :, 2) = CgaCal(q(4, 2, i),q(5, 2, i),q(6, 2, i));
    nj_1=A(:,:,1)*nj1_1;
    nj_2=A(:,:,2)*nj2_2;
    rj(:, 1) = [0; 0; 0];
    rj(:, 2) = q(1:3, 1, i) + A(:, :, 1) * r1_2(:, 1);
    bj(:, 1) = [Vec2Mat(rj(:, 1))*nj_1; nj_1];
    bj(:, 2) = [Vec2Mat(rj(:, 2))*nj_2; nj_2];
    Rd = [bj(:, 1), zeros(6, 1); zeros(6, 1), bj(:, 2)];

    
    temp = T * Rd * dz(:, i);
    dZ(:, 1, i) = temp(1:6);
    dZ(:, 2, i) = temp(7:12);

    wi(:,1)=dZ(4:6, 1, i);
    wi(:,2)=dZ(4:6, 2, i);
    wimat(:,:,1)=Vec2Mat(wi(:,1));
    wimat(:,:,2)=Vec2Mat(wi(:,2));
    dj(:, 1, i) = [(-wimat(:,:,1)*wimat(:,:,1)) * rj(:, 1); wimat(:,:,1)*nj_1*dz(1,i)];
    dj(:, 2, i) = [(wimat(:,:,1)*wimat(:,:,1)-wimat(:,:,2)*wimat(:,:,2)) * rj(:, 2); wimat(:,:,2)*nj_2*dz(2,i)];
    ei(:, 1, i) = [wimat(:,:,1) ^ 2 * q(1:3, 1, i); zeros(3,1)];
    ei(:, 2, i) = [wimat(:,:,2) ^ 2 * q(1:3, 2, i); zeros(3,1)];

    Di(:, :, 1) = [I3, -Vec2Mat(q(1:3, 1, i)); zeros(3), I3];
    Di(:, :, 2) = [I3, -Vec2Mat(q(1:3, 2, i)); zeros(3), I3];
    D = [Di(:, :, 1), zero6; zero6, Di(:, :, 2)];

    Mbar = Rd' * T' * D' * [M(:, :, 1), zero6; zero6, M(:, :, 2)] * D * T * Rd;
    Qebar = Rd' * T' * D' * Qe;
    Qvbar = Rd' * T' * D' * (- [M(:, :, 1), zero6; zero6, M(:, :, 2)] * ([ei(:, 1, i); ei(:, 2, i)] + D * T * [dj(:, 1, i); dj(:, 2, i)]));

    ddz(:, i) = (Mbar ^ -1) * (Qebar + Qvbar);
    dz(:, i + 1) = dz(:, i) + h * ddz(:, i);
    z(:, i + 1) = z(:, i) + h * dz(:, i);

    q(:, 1, i+1) = [-CgaCal(0,0,z(1,i+1)) * r1_1;0;0;z(1,i+1)];
    q(:, 2, i+1) = q(:, 1, i+1) + [(CgaCal(0,0,z(1,i+1)) * r1_2 - CgaCal(0,0,q(6,1,i+1)+z(2,i+1)) * r2_2);0;0;z(2,i+1)];

end

%% figure
figure(1)
theta14fig = reshape(q(6, 1, :), 1, nstep);
theta24fig = reshape(q(6, 2, :), 1, nstep);
hold on
grid on
plot(90 + theta14fig * 180 / pi, '-g');
plot(90 + theta24fig * 180 / pi, '-r');
legend('theta1', 'theta2')
hold off
figure(2)
load('Pendulum.mat');
plot((pi / 2 + theta24fig) * 180 / pi - MEA_ANGLE_1Q);
