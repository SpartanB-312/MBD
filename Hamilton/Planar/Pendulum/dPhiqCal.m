function dPhiq = dPhiqCal(q,dq)
    dPhiq = [0 0 sin(q(3))*dq(3);
            0 0 -cos(q(3))*dq(3)];
end