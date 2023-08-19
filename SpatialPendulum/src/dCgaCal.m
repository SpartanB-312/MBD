function dCGA = dCgaCal(rll,pth,yaw,q)
CGA = [cos(pth)*cos(yaw),     -sin(yaw)*cos(rll)+cos(yaw)*sin(pth)*sin(rll),     sin(yaw)*sin(rll)+cos(yaw)*sin(pth)*cos(rll);
        cos(pth)*sin(yaw),     cos(yaw)*cos(rll)+sin(yaw)*sin(pth)*sin(rll),     -cos(yaw)*sin(rll)+sin(yaw)*sin(pth)*cos(rll);
        -sin(pth) ,     cos(pth)*sin(rll),    cos(pth)*cos(rll)];
%rll
dCGArll=[0,sin(yaw)*sin(rll)+cos(yaw)*sin(pth)*cos(rll),sin(yaw)*cos(rll)-cos(yaw)*sin(pth)*sin(rll);
         0,-cos(yaw)*sin(rll)+sin(yaw)*sin(pth)*cos(rll),-cos(yaw)*cos(rll)-sin(yaw)*sin(pth)*sin(rll);
         0,cos(pth)*cos(rll),-cos(pth)*sin(rll)]*q;
%pth
dCGApth=[-sin(pth)*cos(yaw),     cos(yaw)*cos(pth)*sin(rll),     cos(yaw)*cos(pth)*cos(rll);
        -sin(pth)*sin(yaw),     sin(yaw)*cos(pth)*sin(rll),     sin(yaw)*cos(pth)*cos(rll);
        -cos(pth) ,     -sin(pth)*sin(rll),    -sin(pth)*cos(rll)]*q;
%yaw
dCGAyaw=[-cos(pth)*sin(yaw),     -cos(yaw)*cos(rll)+-sin(yaw)*sin(pth)*sin(rll),     cos(yaw)*sin(rll)-sin(yaw)*sin(pth)*cos(rll);
        cos(pth)*cos(yaw),     -sin(yaw)*cos(rll)+cos(yaw)*sin(pth)*sin(rll),     sin(yaw)*sin(rll)+cos(yaw)*sin(pth)*cos(rll);
        0 ,     0,    0]*q;
%
dCGA=[dCGArll,dCGApth,dCGAyaw];