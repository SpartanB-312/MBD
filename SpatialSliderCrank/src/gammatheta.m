function P = gammatheta(Ai,Aj,wi,wj,hi)
    wiMat=Vec2Mat(wi);
    wjMat=Vec2Mat(wj);
    P=-hi'*(Ai'*Aj*wjMat-wiMat*Ai'*Aj)*wj;
end