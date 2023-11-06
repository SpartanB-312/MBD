function K = KCal(e,lam,mu,A,l)
    K1 = K1Cal(lam,mu,A,l);
    K2 = K2Cal(e,lam,mu,A,l);
    K=K1+K2;
end