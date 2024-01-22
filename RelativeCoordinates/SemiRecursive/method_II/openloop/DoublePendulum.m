%4.9
clear;
m1 = 1; M(:, :, 1) = diag([1, 1, 1/4]);
m2 = 1; M(:, :, 2) = diag([1, 1, 1/4]);
g = -9.8;
m1 = 1; M(:, :, 1) = diag([1, 1, 0]);
m2 = 1; M(:, :, 2) = diag([1, 1, 0]);

I2 = eye(2); I3 = eye(3); zero3 = zeros(3);
I2tilde = [0 -1; 1 0];

T = [I3, zero3; I3, I3];

r1_1 = [-1; 0];
r1_2 = [0; 0];
r2_2 = [-1; 0];

Z(:, 1, 1) = [0; 0; -pi / 3];
Z(:, 2, 1) = [0; 0; -pi / 3];
q(:, 1, 1) = [0.25; -sqrt(3) / 4; -pi / 3];
q(:, 2, 1) = [1.5/2; -1.5 * sqrt(3) / 2; -pi / 3];
%q(:,2,1)=[0.5;-0.5*sqrt(3)/2-0.25;-pi/2];
q(:, 1, 1) = [0.5; -sqrt(3) / 2; -pi / 3];
q(:, 2, 1) = [1; -sqrt(3); -pi / 3];
A(:, :, 1) = CgaCal(q(3, 1, 1));
A(:, :, 2) = CgaCal(q(3, 2, 1));
rj(:, 1) = [0; 0];
rj(:, 2) = q(1:2, 1, 1) + A(:, :, 1) * r1_2(:, 1);
bj(:, 1) = [-I2tilde * rj(:, 1); 1];
bj(:, 2) = [-I2tilde * rj(:, 2); 1];
Rd = [bj(:, 1), zeros(3, 1); zeros(3, 1), bj(:, 2)];

%Z(:, 1, 1) = [CgaCal(z(1,1)) * r1_1;z(1,1)];
%Z(:, 2, 1) = Z(:, 1, 1) + [(CgaCal(z(1,1)) * r1_2 - CgaCal(Z(3,1,1)+z(2,1)) * r2_2);z(2,1)];

dz(:, 1) = [0; 0];
dj(:, 1, 1) = [(dz(1, 1) ^ 2) * rj(:, 1); 0];
dj(:, 2, 1) = [(dz(2, 1) ^ 2 - dz(1, 1) ^ 2) * rj(:, 2); 0];
temp = T * Rd * dz(:, 1);
dZ(:, 1, 1) = temp(1:3);
dZ(:, 2, 1) = temp(4:6);
%temp = T * (Rd * ddz(:, 1) + [dj(:, 1, 1); dj(:, 2, 1)]);
%ddZ(:, 1, 1) = temp(1:3);
%ddZ(:, 2, 1) = temp(4:6);


ei(:, 1, 1) = [-dZ(3, 1, 1) ^ 2 * q(1:2, 1, 1); 0];
ei(:, 2, 1) = [-dZ(3, 2, 1) ^ 2 * q(1:2, 2, 1); 0];

Di(:, :, 1) = [I2, I2tilde * q(1:2, 1, 1); zeros(1, 2), 1];
Di(:, :, 2) = [I2, I2tilde * q(1:2, 2, 1); zeros(1, 2), 1];
D = [Di(:, :, 1), zero3; zero3, Di(:, :, 2)];
dq(:, 1, 1) = Di(:, :, 1) * dZ(:, 1, 1);
dq(:, 2, 1) = Di(:, :, 1) * dZ(:, 2, 1);

z(:, 1) = [-pi / 3; 0];

Qe(1:3, 1) = [0; m1 * g; 0];
Qe(4:6, 1) = [0; m2 * g; 0];

Time = 5; h = 0.001; nstep = Time / h;

for i = 1:nstep - 1

    A(:, :, 1) = CgaCal(q(3, 1, i));
    A(:, :, 2) = CgaCal(q(3, 2, i));
    rj(:, 1) = [0; 0];
    rj(:, 2) = q(1:2, 1, i) + A(:, :, 1) * r1_2;
    bj(:, 1) = [-I2tilde * rj(:, 1); 1];
    bj(:, 2) = [-I2tilde * rj(:, 2); 1];
    Rd = [bj(:, 1), zeros(3, 1); zeros(3, 1), bj(:, 2)];

    
    temp = T * Rd * dz(:, i);
    dZ(:, 1, i) = temp(1:3);
    dZ(:, 2, i) = temp(4:6);
    dj(:, 1, i) = [(dZ(3,1, i) ^ 2) * rj(:, 1); 0];
    dj(:, 2, i) = [(dZ(3,2, i) ^ 2 - dZ(3,1, i) ^ 2) * rj(:, 2); 0];

    ei(:, 1, i) = [-dZ(3, 1, i) ^ 2 * q(1:2, 1, i); 0];
    ei(:, 2, i) = [-dZ(3, 2, i) ^ 2 * q(1:2, 2, i); 0];

    Di(:, :, 1) = [I2, I2tilde * q(1:2, 1, i); zeros(1, 2), 1];
    Di(:, :, 2) = [I2, I2tilde * q(1:2, 2, i); zeros(1, 2), 1];
    D = [Di(:, :, 1), zero3; zero3, Di(:, :, 2)];

    Mbar = Rd' * T' * D' * [M(:, :, 1), zeros(3); zeros(3), M(:, :, 2)] * D * T * Rd;
    Qebar = Rd' * T' * D' * Qe;
    Qvbar = Rd' * T' * D' * (- [M(:, :, 1), zeros(3); zeros(3), M(:, :, 2)] * ([ei(:, 1, i); ei(:, 2, i)] + D * T * [dj(:, 1, i); dj(:, 2, i)]));

    ddz(:, i) = (Mbar ^ -1) * (Qebar + Qvbar);
    dz(:, i + 1) = dz(:, i) + h * ddz(:, i);
    z(:, i + 1) = z(:, i) + h * dz(:, i);

    q(:, 1, i+1) = [-CgaCal(z(1,i+1)) * r1_1;z(1,i+1)];
    q(:, 2, i+1) = q(:, 1, i+1) + [(CgaCal(z(1,i+1)) * r1_2 - CgaCal(q(3,1,i+1)+z(2,i+1)) * r2_2);z(2,i+1)];

end

%% figure
figure(1)
theta14fig = reshape(q(3, 1, :), 1, nstep);
theta24fig = reshape(q(3, 2, :), 1, nstep);
hold on
grid on
plot(90 + theta14fig * 180 / pi, '-g');
plot(90 + theta24fig * 180 / pi, '-r');
legend('theta1', 'theta2')
hold off
figure(2)
load('DoublePendulum.mat');
plot((pi / 2 + theta14fig) * 180 / pi - MEA_ANGLE_1Q);
