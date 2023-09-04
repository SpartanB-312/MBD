function D = AngleRateDCal(rll,pth)
    D = [1,  sin(rll)*tan(pth),    cos(rll)*tan(pth)     
        0,  cos(rll),              -sin(rll)           ;
        0,  sin(rll)/cos(pth),    cos(rll)/cos(pth)];
end