function M = MassCal(rho,l)
    syms x;
    I2=[1,0;0,1];
    S1=1-3*x^2+2*x^3;
    S2=l*(x-2*x^2+x^3);
    S3=3*x^2-2*x^3;
    S4=-l*(x^2-x^3);
    S=[I2*S1,I2*S2,I2*S3,I2*S4];
    f=S'*S;
    M=rho*eval(int(f,x,0,1));
end