function constr = constraint_function_fractional(parS, h, C, g0, galpha, g1)
    
    gam = parS.gam;
    
%     x  = g0 * C(1:h);
    x1 = g1 * C(1:h);
%     x2 = g2 * C(1:h);
    
    y  = g0 * C(h+1:2*h);
    yalpha = galpha * C(h+1:2*h);
%     y1 = g1 * C(h+1:2*h);
    
    constr = yalpha + gam*y - x1;

end

