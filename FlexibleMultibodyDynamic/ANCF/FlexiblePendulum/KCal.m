function K = KCal(e,lam,mu,h,l,w)
    K1 = K1Cal(lam,mu,h,l,w);
    K2 = K2Cal(e,lam,mu,h,l,w);
    K=K1+K2;
end