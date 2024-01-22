function P = gammad1(ai,aj,Ai,Aj,wi,wj)
    wiMat=Vec2Mat(wi);
    wjMat=Vec2Mat(wj);
    aiMat=Vec2Mat(ai);
    ajMat=Vec2Mat(aj);
    P=-aj'*(Aj'*Ai*wiMat*wiMat+wjMat*wjMat*Aj'*Ai)*ai+2*wj'*ajMat*Aj'*Ai*aiMat*wi;
end