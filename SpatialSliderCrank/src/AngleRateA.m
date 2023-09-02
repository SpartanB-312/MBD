function [DYAW, DPTH, DRLL] =  AngleRateA(p,q,r,rll,pth)
temp = [0,  sin(rll)/cos(pth),    cos(rll)/cos(pth)    ;
        0,  cos(rll),              -sin(rll)           ;
        1,  sin(rll)*tan(pth),    cos(rll)*tan(pth)     ]*[p;q;r];
DYAW = temp(1);
DPTH = temp(2);
DRLL = temp(3);