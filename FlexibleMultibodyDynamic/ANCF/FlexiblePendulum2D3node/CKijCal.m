function CKij = CKijCal(lam,mu,h,l,w)
    CKijtemp=zeros(12,12,12,12);
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
    Sxx=(Sx'*Sx);
    Sxy=(Sx'*Sy);
    Syx=(Sy'*Sx);
    Syy=(Sy'*Sy);
    for i = 1:12
        for j = 1:12
            f1xy=Sxx(i,:)'*Sxx(j,:)+Syy(i,:)'*Syy(j,:);
            f1y=int(f1xy,x,0,l);
            intf1=int(f1y,y,-w/2,w/2);

            f2xy=Sxx(i,:)'*Syy(j,:)+Syy(i,:)'*Sxx(j,:);
            f2y=int(f2xy,x,0,l);
            intf2=int(f2y,y,-w/2,w/2);

            f3xy=Sxy(i,:)'*Syx(j,:)+Syx(i,:)'*Sxy(j,:);
            f3y=int(f3xy,x,0,l);
            intf3=int(f3y,y,-w/2,w/2);

            CKijtemp(:,:,i,j)=(lam+2*mu)/2*intf1+lam/2*intf2+mu*intf3;
        end
    end
    
    CKij=CKijtemp*h;
end