function K1 = K1Cal(lam,mu,A,l)
    syms x;
    I2=[1,0;0,1];
    S1=1-3*x^2+2*x^3;
    S2=l*(x-2*x^2+x^3);
    S3=3*x^2-2*x^3;
    S4=-l*(x^2-x^3);
    S=[I2*S1,I2*S2,I2*S3,I2*S4];
    Sx=diff(S,x);
    f=Sx'*Sx*A;
    intf=int(f,x,0,1);
    K1=-(lam+mu)*eval(intf);
end