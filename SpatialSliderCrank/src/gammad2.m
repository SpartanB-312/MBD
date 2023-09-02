function P = gammad2(ai,Ai,Aj,wi,wj,vi,vj,si,sj,dij)
    wiMat=Vec2Mat(wi);
    wjMat=Vec2Mat(wj);
    aiMat=Vec2Mat(ai);
    P=2*wi'*aiMat*Ai'*(vj-vi) + 2*sj'*wjMat*Aj'*Ai*wiMat*ai - ...
        si'*wiMat*wiMat*ai - sj'*wjMat*wjMat*Aj'*Ai*ai - dij'*Ai*wiMat*wiMat*ai;
end