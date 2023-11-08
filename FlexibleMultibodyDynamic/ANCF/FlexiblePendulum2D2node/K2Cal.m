function K2 = K2Cal(e,lam,mu,h,l,w)
    K2temp=zeros(12);
    CKij = CKijCal(lam,mu,h,l,w);
    for i = 1:12
        for j = 1:12
            K2temp(i,j)=e'*CKij(:,:,i,j)*e;
        end
    end
    K2=K2temp;
end