function K2 = K2Cal(e,lam,mu,A,l)
    K2temp=zeros(8);
    CKij = CKijCal(lam,mu,A,l);
    for i = 1:8
        for j = 1:8
            K2temp(i,j)=e'*CKij(:,:,i,j)*e;
        end
    end
    K2=K2temp;
end