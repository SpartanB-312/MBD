function Phiq = PhiqCal(q)
    Phiq = [1 0 -cos(q(3));
            0 1 -sin(q(3))];
end