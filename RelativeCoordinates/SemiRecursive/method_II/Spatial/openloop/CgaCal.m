function CGA = CgaCal(rll,pth,yaw)
CGA = [cos(pth)*cos(yaw),     -sin(yaw)*cos(rll)+cos(yaw)*sin(pth)*sin(rll),     sin(yaw)*sin(rll)+cos(yaw)*sin(pth)*cos(rll);
    cos(pth)*sin(yaw),     cos(yaw)*cos(rll)+sin(yaw)*sin(pth)*sin(rll),     -cos(yaw)*sin(rll)+sin(yaw)*sin(pth)*cos(rll);
   -sin(pth) ,     cos(pth)*sin(rll),    cos(pth)*cos(rll)];