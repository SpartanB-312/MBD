%4.9
clear;
g = -9.8;
m1 = 0; M(:, :, 1) = diag([m1, m1, m1, 0, 0, 0]);
m2 = 0; M(:, :, 2) = diag([m2, m2, m2, 0, 0, 0]);
m3 = 0; M(:, :, 3) = diag([m3, m3, m3, 0, 0, 0]);
m4 = 0; M(:, :, 4) = diag([m4, m4, m4, 0, 0, 0]);
m5 = 0; M(:, :, 5) = diag([m5, m5, m5, 0, 0, 0]);
m6 = 1; J6=1; M(:, :, 6) = diag([m6, m6, m6, J6, J6, J6]);

I2 = eye(2); I3 = eye(3); I6 = eye(6);
zero3 = zeros(3); zero6 = zeros(6);
I2tilde = [0 -1; 1 0];

T = [I6, zero6, zero6, zero6, zero6, zero6;
     I6, I6, zero6, zero6, zero6, zero6;
     I6, I6, I6, zero6, zero6, zero6;
     I6, I6, I6, I6, zero6, zero6;
     I6, I6, I6, I6, I6, zero6;
     I6, I6, I6, I6, I6, I6];

v_L=0;%虚拟长度
r1_1 = [0; 0; 0];r1_2 = [0; 0; 0];
r2_2 = [0; 0; 0];r2_3 = [0; 0; 0];
r3_3 = [0; 0; 0];r3_4 = [0; 0; 0];
r4_4 = [0; 0; 0];r4_5 = [0; 0; 0];
r5_5 = [0; 0; 0];r5_6 = [0; 0; 0];
r6_6 = [0; 0; 0];
%平移轴
nj1_1 = [1; 0; 0];nj1_2 = [0; 1; 0];
nj2_2 = [0; 1; 0];nj2_3 = [0; 0; 1];
nj3_3 = [0; 0; 1];nj3_4 = [0; 0; 1];
%旋转轴zyx
nj4_4 = [0; 0; 1];nj4_5 = [0; 1; 0];
nj5_5 = [0; 1; 0];nj5_6 = [1; 0; 0];
nj6_6 = [1; 0; 0];

for j = 1:6
    q(:,j,1)=zeros(6,1);
    A(:, :, j) = CgaCal(q(4, j, 1),q(5, j, 1),q(6, j, 1));
end
z(:, 1) = [0; 0; 0; 0; 0; 0];
q(1,6,1)=0;

nj_1=A(:,:,1)*nj1_1;
nj_2=A(:,:,2)*nj2_2;
nj_3=A(:,:,3)*nj3_3;
nj_4=A(:,:,4)*nj4_4;
nj_5=A(:,:,5)*nj5_5;
nj_6=A(:,:,6)*nj6_6;

rj(:, 1) = [0; 0; 0];
for j = 2:6
    rj(:, j) = q(1:3, j-1, 1) + A(:, :, j-1) * r1_2;
end

%平移
bj(:, 1) = [nj_1; zeros(3,1)];
bj(:, 2) = [nj_2; zeros(3,1)];
bj(:, 3) = [nj_3; zeros(3,1)];
%旋转
bj(:, 4) = [Vec2Mat(rj(:, 4))*nj_4; nj_4];
bj(:, 5) = [Vec2Mat(rj(:, 5))*nj_5; nj_5];
bj(:, 6) = [Vec2Mat(rj(:, 6))*nj_6; nj_6];
Rd = [bj(:, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1);
      zeros(6, 1), bj(:, 2), zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1);
      zeros(6, 1), zeros(6, 1), bj(:, 3), zeros(6, 1), zeros(6, 1), zeros(6, 1);
      zeros(6, 1), zeros(6, 1), zeros(6, 1), bj(:, 4), zeros(6, 1), zeros(6, 1);
      zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1), bj(:, 5), zeros(6, 1);
      zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1), bj(:, 6)];


dz(:, 1) = [0; 0; 0; 0; 0; 0];
temp = T * Rd * dz(:, 1);
for j = 1:6
    dZ(:,j,1)=temp(1+6*(j-1):6+6*(j-1));
    wi(:,j)=dZ(4:6, j, 1);
    wimat(:,:,j)=Vec2Mat(wi(:,j));
end

%
dj(:, 1, 1) = [wimat(:,:,1)*nj_1*dz(1,1);zeros(3,1)];
dj(:, 2, 1) = [wimat(:,:,2)*nj_2*dz(2,1);zeros(3,1)];
dj(:, 3, 1) = [wimat(:,:,3)*nj_3*dz(3,1);zeros(3,1)];
%
dj(:, 4, 1) = [(wimat(:,:,3)*wimat(:,:,3)-wimat(:,:,4)*wimat(:,:,4)) * rj(:, 4); wimat(:,:,4)*nj_4*dz(4,1)];
dj(:, 5, 1) = [(wimat(:,:,4)*wimat(:,:,4)-wimat(:,:,5)*wimat(:,:,5)) * rj(:, 5); wimat(:,:,5)*nj_5*dz(5,1)];
dj(:, 6, 1) = [(wimat(:,:,5)*wimat(:,:,5)-wimat(:,:,6)*wimat(:,:,6)) * rj(:, 6); wimat(:,:,6)*nj_6*dz(6,1)];

for j = 1:6
    ei(:, j, 1) = [wimat(:,:,j) ^ 2 * q(1:3, j, 1); zeros(3,1)];
    Di(:, :, j) = [I3, -Vec2Mat(q(1:3, j, 1)); zeros(3), I3];
    dq(:, j, 1) = Di(:, :, j) * dZ(:, j, 1);
end

D = [Di(:, :, 1), zero6, zero6, zero6, zero6, zero6;
     zero6, Di(:, :, 2), zero6, zero6, zero6, zero6;
     zero6, zero6, Di(:, :, 3), zero6, zero6, zero6;
     zero6, zero6, zero6, Di(:, :, 4), zero6, zero6;
     zero6, zero6, zero6, zero6, Di(:, :, 5), zero6;
     zero6, zero6, zero6, zero6, zero6, Di(:, :, 6)];

Qe(:, 1) = zeros(36,1);
Qe(31) = 1;Qe(32)=0;Qe(33)=0;
Qe(34) = 1;Qe(35)=0;Qe(36)=0;

Time = 2; h = 0.001; nstep = Time / h;

for i = 1:nstep - 1

    for j = 1:6
        A(:, :, j) = CgaCal(q(4, j, i),q(5, j, i),q(6, j, i));
    end

    nj_1=A(:,:,1)*nj1_1;
    nj_2=A(:,:,2)*nj2_2;
    nj_3=A(:,:,3)*nj3_3;
    nj_4=A(:,:,4)*nj4_4;
    nj_5=A(:,:,5)*nj5_5;
    nj_6=A(:,:,6)*nj6_6;

    rj(:, 1) = [0; 0; 0];
    for j = 2:6
        rj(:, j) = q(1:3, j-1, i) + A(:, :, j-1) * r1_2;
    end
    %平移
    bj(:, 1) = [nj_1; zeros(3,1)];
    bj(:, 2) = [nj_2; zeros(3,1)];
    bj(:, 3) = [nj_3; zeros(3,1)];
    %旋转
    bj(:, 4) = [Vec2Mat(rj(:, 4))*nj_4; nj_4];
    bj(:, 5) = [Vec2Mat(rj(:, 5))*nj_5; nj_5];
    bj(:, 6) = [Vec2Mat(rj(:, 6))*nj_6; nj_6];

    Rd = [bj(:, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1);
          zeros(6, 1), bj(:, 2), zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1);
          zeros(6, 1), zeros(6, 1), bj(:, 3), zeros(6, 1), zeros(6, 1), zeros(6, 1);
          zeros(6, 1), zeros(6, 1), zeros(6, 1), bj(:, 4), zeros(6, 1), zeros(6, 1);
          zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1), bj(:, 5), zeros(6, 1);
          zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1), bj(:, 6)];

    
    temp = T * Rd * dz(:, i);
    for j = 1:6
        dZ(:,j,i)=temp(1+6*(j-1):6+6*(j-1));
        wi(:,j)=dZ(4:6, j, i);
        wimat(:,:,j)=Vec2Mat(wi(:,j));
    end

    %
    dj(:, 1, i) = [wimat(:,:,1)*nj_1*dz(1,i);zeros(3,1)];
    dj(:, 2, i) = [wimat(:,:,2)*nj_2*dz(2,i);zeros(3,1)];
    dj(:, 3, i) = [wimat(:,:,3)*nj_3*dz(3,i);zeros(3,1)];
    %
    dj(:, 4, i) = [(wimat(:,:,3)*wimat(:,:,3)-wimat(:,:,4)*wimat(:,:,4)) * rj(:, 4); wimat(:,:,4)*nj_4*dz(4,i)];
    dj(:, 5, i) = [(wimat(:,:,4)*wimat(:,:,4)-wimat(:,:,5)*wimat(:,:,5)) * rj(:, 5); wimat(:,:,5)*nj_5*dz(5,i)];
    dj(:, 6, i) = [(wimat(:,:,5)*wimat(:,:,5)-wimat(:,:,6)*wimat(:,:,6)) * rj(:, 6); wimat(:,:,6)*nj_6*dz(6,i)];

    for j = 1:6
        ei(:, j, i) = [wimat(:,:,j) ^ 2 * q(1:3, j, i); zeros(3,1)];
        Di(:, :, j) = [I3, -Vec2Mat(q(1:3, j, i)); zeros(3), I3];
        dq(:, j, i) = Di(:, :, j) * dZ(:, j, i);
    end
    D = [Di(:, :, 1), zero6, zero6, zero6, zero6, zero6;
         zero6, Di(:, :, 2), zero6, zero6, zero6, zero6;
         zero6, zero6, Di(:, :, 3), zero6, zero6, zero6;
         zero6, zero6, zero6, Di(:, :, 4), zero6, zero6;
         zero6, zero6, zero6, zero6, Di(:, :, 5), zero6;
         zero6, zero6, zero6, zero6, zero6, Di(:, :, 6)];

    Mdiag = [M(:, :, 1), zero6, zero6, zero6, zero6, zero6;
             zero6, M(:, :, 2), zero6, zero6, zero6, zero6;
             zero6, zero6, M(:, :, 3), zero6, zero6, zero6;
             zero6, zero6, zero6, M(:, :, 4), zero6, zero6;
             zero6, zero6, zero6, zero6, M(:, :, 5), zero6;
             zero6, zero6, zero6, zero6, zero6, M(:, :, 6)];

    ei4Qvbar=ei(:,1,i);
    dj4Qvbar=dj(:,1,i);
    for j = 2:6
        ei4Qvbar=[ei4Qvbar;ei(:, j, i)];
        dj4Qvbar=[dj4Qvbar;dj(:, j, i)];
    end

    Mbar = Rd' * T' * D' * Mdiag * D * T * Rd;
    Qebar = Rd' * T' * D' * Qe;
    Qvbar = Rd' * T' * D' * (- Mdiag * (ei4Qvbar + D * T * dj4Qvbar));

    ddz(:, i) = Mbar \ (Qebar + Qvbar);
    dz(:, i + 1) = dz(:, i) + h * ddz(:, i);
    z(:, i + 1) = z(:, i) + h * dz(:, i);

    q(:, 1, i+1) = [z(1,i+1);0;0;0;0;0];
    q(:, 2, i+1) = q(:, 1, i+1) + [0;z(2,i+1);0;0;0;0];
    q(:, 3, i+1) = q(:, 2, i+1) + [0;0;z(3,i+1);0;0;0];
    q(:, 4, i+1) = q(:, 3, i+1) + [0;0;0;0;0;z(4,i+1)];
    q(:, 5, i+1) = q(:, 4, i+1) + [0;0;0;0;z(5,i+1);0];
    q(:, 6, i+1) = q(:, 5, i+1) + [0;0;0;z(6,i+1);0;0];

end

%% figure
figure(1)
x4fig = reshape(q(1, 6, :), 1, nstep);
y4fig = reshape(q(2, 6, :), 1, nstep);
z4fig = reshape(q(3, 6, :), 1, nstep);
hold on
grid on
plot(x4fig, '-g');
plot(y4fig, '-r');
plot(z4fig, '-b');
legend('x', 'y','z')
hold off
figure(2)
xtheta4fig = reshape(q(4, 6, :), 1, nstep);
grid on
plot(xtheta4fig, '-g');