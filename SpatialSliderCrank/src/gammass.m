function P = gammass(Ai,Aj,wi,wj,vi,vj,si,sj,dij)
    wiMat=Vec2Mat(wi);
    wjMat=Vec2Mat(wj);
    siMat=Vec2Mat(si);
    sjMat=Vec2Mat(sj);
    temp1=(vj+Aj*wjMat*sj-vi-Ai*wiMat*si)';
    temp2=(vi-vj)-Ai*siMat*wi+Aj*sjMat*wj;
    temp3=dij'*(Ai*wiMat*siMat*wi-Aj*wjMat*sjMat*wj);
    P=2*temp1*temp2-2*temp3;
end