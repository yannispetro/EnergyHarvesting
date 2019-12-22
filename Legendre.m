function [g0,galpha,g1,g2] = Legendre(parW,Calpha)
    
    ti = parW.ti;
    tf = parW.tf;
    Fs = parW.Fs;
    h  = parW.h; 
    
    
    % - Time discretization and formation of the matrices
    dt = 1/Fs;                          % Sampling period (Original - sec)
    nt  = floor((tf - ti)/dt) + 1;      % Length of signal
    tt = (ti + (0:nt-1)*dt).';          % Time vector
    
    % - Construction shifted Legendre polynomials -
    syms t
    P = sym('P', [1 h]);
    P(1) = 1;
    P(2) = (2 * t - ti - tf) / (tf - ti);

    for kk = 3:h
        p = kk - 2;
        P(kk) = expand((2 * p + 1) / (p + 1) * (2 * t - ti - tf)...
            / (tf - ti) * P(kk - 1) - p / (p + 1) * P(kk - 2));
    end
    clear kk p
    
    syms time
    % - g is the function that each of its elements multiplies one c_i 
%     g = ((t - ti) ^ ord) * ((t - tf) ^ ord) * P;
    g = P;
    if Calpha < 1
        alphadgdt = (1/gamma(1-Calpha))*int(diff(g,t)/(time-t)^Calpha,t,ti,time); %-- Change due to fractional
        alphadgdt = subs(alphadgdt, time, t);
    end
    dgdt = diff(g, t);
    d2gdt2 = diff(g, t, 2);

    % convert g and derivatives to matrices
    g0 = double(subs(g, tt));
    g1 = double(subs(dgdt, tt));
    g2 = double(subs(d2gdt2, tt));
    if Calpha < 1
        galpha = double(subs(alphadgdt, tt)); %-- Change due to fractional
    else
        galpha = g1;
    end
end

