function df = Derivation(f0,f1,dq,sizedq)
    %f0=f(pars0);
    %pars1=pars0+[zeros(nf),dq',zeros(nb)];
    [sizef0,~]=size(f0);
    df=zeros(sizef0,sizedq);

    for i=1:sizedq
        %pars1=pars0+[zeros(1,nf+i-1),dq,zeros(nb+sizedq-i)];
        %f1=f(pars1);
        df(:,i)=(f1(:,i)-f0)/dq;
    end

end