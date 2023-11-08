function K1 = K1Cal(lam,mu,h,l,w)
    syms x y;
    I2=[1,0;0,1];
    S1=1-3*(x/l)+2*(x/l)^2;
    S2=y*(1-3*(x/l)+2*(x/l)^2);
    S3=4*(x/l-(x/l)^2);
    S4=y*(4*(x/l)-4*(x/l)^2);
    S5=-x/l+2*(x/l)^2;
    S6=y*(-x/l+2*(x/l)^2);
    S=[I2*S1,I2*S2,I2*S3,I2*S4,I2*S5,I2*S6];
    Sx=diff(S,x);
    Sy=diff(S,y);
    fxy=(Sx'*Sx+Sy'*Sy);
    fy=int(fxy,x,0,l);
    intf=int(fy,y,-w/2,w/2);
    K1=-(lam+mu)*eval(intf)*h;
end