function Qg = QgCal(rho,g,h,l,w)
    syms x y;
    I2=[1,0;0,1];
    S1 = 1 - 3*(x/l)^2 + 2*(x/l)^3;
    S2 = l*(x/l - 2*(x/l)^2 + (x/l)^3);
    S3 = l*(y/l - (x*y)/(l^2));
    S4 = 3*(x/l)^2 - 2*(x/l)^3;
    S5 = l*(-(x/l)^2 + (x/l)^3);
    S6 = l*((x*y)/(l^2));
    S=[I2*S1,I2*S2,I2*S3,I2*S4,I2*S5,I2*S6];
    fxy=S(2,:)';
    fy=int(fxy,x,0,l);
    intf=int(fy,y,-w/2,w/2);
    Qg=rho*g*eval(intf)*h;
end