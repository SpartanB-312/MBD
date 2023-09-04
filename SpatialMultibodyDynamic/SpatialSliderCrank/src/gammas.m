function P = gammas(Ai,Aj,wi,wj,si,sj)
    wiMat=Vec2Mat(wi);
    wjMat=Vec2Mat(wj);
    P=Ai*wiMat*wiMat*si-Aj*wjMat*wjMat*sj;
end