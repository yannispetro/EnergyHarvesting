function J = stochastic_action(parS, h, C, g0, g1, g2, tt)

    zeta   = parS.zeta;
    lambda = parS.lambda;
    delta  = parS.delta;
    kappa  = parS.kappa;
%     gam    = parS.gam;
    S0     = parS.S0;

    x  = g0 * C(1:h);
    x1 = g1 * C(1:h);
    x2 = g2 * C(1:h);
    
    y  = g0 * C(h+1:2*h);
%     y1 = g1 * C(h+1:2*h);

    sys = (x2 + 2*zeta*x1 + x + lambda*x.^2 + delta*x.^3 + kappa^2*y);
    
    L = 1/(4*pi*S0)*( sys ).^2;

    J = trapz(tt,L);
end
