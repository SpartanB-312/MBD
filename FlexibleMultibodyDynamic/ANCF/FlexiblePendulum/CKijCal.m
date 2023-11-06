function CKij = CKijCal(lam,mu,A,l)
    CKijtemp=zeros(8,8,8,8);
    syms x;
    I2=[1,0;0,1];
    S1=1-3*x^2+2*x^3;
    S2=l*(x-2*x^2+x^3);
    S3=3*x^2-2*x^3;
    S4=-l*(x^2-x^3);
    S=[I2*S1,I2*S2,I2*S3,I2*S4];
    Sx=diff(S,x);
    S11=(Sx'*Sx);
    for i = 1:8
        for j = 1:8
            f1=S11(i,:)'*S11(j,:)*A;
            intf1=int(f1,x,0,1);

            CKijtemp(:,:,i,j)=(lam+2*mu)/2*intf1;
        end
    end
    
    CKij=CKijtemp;
end