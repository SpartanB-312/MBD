function [rll,pth,yaw] = EulerAngleCal(xyzLocal,xyzGlobal)
    %B=RA
    %xyzGlobal=R*xyzLocal
    C=cross(xyzLocal,xyzGlobal);
    theta=acos(xyzLocal'*xyzGlobal/(norm(xyzLocal)*norm(xyzGlobal)));
    R11=C(1)*C(1)*(1-cos(theta))+cos(theta);
    R21=C(1)*C(2)*(1-cos(theta))+C(3)*sin(theta);
    R31=C(1)*C(3)*(1-cos(theta))-C(2)*sin(theta);
    R32=C(2)*C(3)*(1-cos(theta))+C(1)*sin(theta);
    R33=C(3)*C(3)*(1-cos(theta))+cos(theta);

    pth=atan2(-R31,sqrt(R11^2+R21^2));
    yaw=atan2(R21/cos(pth),R11/cos(pth));
    rll=atan2(R32/cos(pth),R33/cos(pth));

end