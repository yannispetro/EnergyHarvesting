function [PDF, Xgrid] = a2_WPI_function_opt(parS, parW, parO )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Calpha = parS.Calpha;

    ti     = parW.ti;
    tf     = parW.tf;
    points = parW.points;
%     domain = parW.domain;
    Fs     = parW.Fs;
    h      = parW.h;
    targetPDF = parW.targetPDF;
    g0 = parW.g0;
    galpha = parW.galpha;
    g1 = parW.g1;
    g2 = parW.g2;

    TolsIP  = parO.TolsIP;
    MaxIter = parO.MaxIter;
    MaxFun  = parO.MaxFun;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Addapt domain
    scal = 1;
    Zq = [parS.gam,parS.delta,parS.zeta];
    if parS.gam > 2.0
        Zq(1) = 2.0;
    end
    if parS.delta > 1.5
        Zq(2) = 1.5;
    end
    load('Model.mat')
    if targetPDF == 1
        XLq = scal*griddatan(Model.data,Model.XL,Zq,'linear');
        XUq = scal*griddatan(Model.data,Model.XU,Zq,'linear');

        if parS.delta == 0
            xlu = max(abs(XLq),abs(XUq));
            XLq2 = -xlu;
            XUq2 = xlu;
        else
            cntr = ( XLq + XUq )/2;
            smtr = (2*sqrt(parS.delta) - parS.lambda)/(2*sqrt(parS.delta));
            XLq2 = XLq - cntr*smtr;
            XUq2 = XUq - cntr*smtr;
            if parS.delta > 1.5
                XLq2 = XLq2*(1-(parS.delta-1.5)*0.35/3.5);
                XUq2 = XUq2*(1-(parS.delta-1.5)*0.35/3.5);
            end
        end
        Xgrid = linspace(XLq2, XUq2, parW.points);
%         Xgrid = linspace(XLq, XUq, parW.points);

    elseif targetPDF == 3
        YLq = scal*griddatan(Model.data,Model.YL,Zq,'linear');
        YUq = scal*griddatan(Model.data,Model.YU,Zq,'linear');
        
        ylu = max(abs(YLq),abs(YUq));
        
        if parS.gam > 2
            ylu = ylu*(1-(parS.gam-2.0)*0.35);
        end
        
        Xgrid = linspace(-ylu, ylu, parW.points);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % - Time discretization and formation of the matrices
    dt = 1/Fs;                          % Sampling period (Original - sec)
    nt  = floor((tf - ti)/dt) + 1;      % Length of signal
    tt = (ti + (0:nt-1)*dt).';          % Time vector
    
    % Nullspace for consideration of BCs
    if targetPDF == 1
        Acon = [g0(1,:),zeros(1,h);g1(1,:),zeros(1,h);zeros(1,h),g0(1,:);g0(end,:),zeros(1,h)];
    elseif targetPDF == 2
        Acon = [g0(1,:),zeros(1,h);g1(1,:),zeros(1,h);zeros(1,h),g0(1,:);g1(end,:),zeros(1,h)];
    elseif targetPDF == 3
        Acon = [g0(1,:),zeros(1,h);g1(1,:),zeros(1,h);zeros(1,h),g0(1,:);zeros(1,h),g0(end,:)];
    end
%     Z = null(Acon);
%     Nvar = size(Z,2);
%     size(Acon)
%     size(Z)
    
    % definition of objective and constraint functions
    func = @(C) stochastic_action(parS, h, C, g0, g1, g2, tt);
    if Calpha < 1
        Constraint_integral = @(C) trapz(tt,(constraint_function_fractional(parS, h, C, g0, galpha, g1)).^2 );
    else
        Constraint_integral = @(C) trapz(tt,(constraint_function_ordinary(parS, h, C, g0, g1)).^2 );
    end
        
    % --- MAIN LOOP ----

    % Create the original grid
    Np = points;

    actions = zeros([1, Np]);
    constraints = zeros([1, Np]);
    exit_flags = zeros([1, Np]);
    C0 = zeros(2*h,1);
    opts = optimoptions('fmincon','Display','none','Algorithm','interior-point',...
            'MaxIterations',MaxIter,...
            'MaxFunctionEvaluations',MaxFun,...
            'OptimalityTolerance',TolsIP(1),...
            'StepTolerance',TolsIP(2),...
            'ConstraintTolerance',TolsIP(3));
    for r = 1:Np
        warning('off')

        bcon = [0;0;0;Xgrid(r)];
        
        
        objfun  = @(C) func(C);
        confun  = @(C) Constraint_integral(C);
        [X,fval,exitflag,output] = fmincon(objfun,C0,[],[],Acon,bcon,[],[],@(x)nonlcon(x,confun),opts);
        C0 = X;
        
        actions(r) = objfun(X);
        constraints(r) = confun(X);
        exit_flags(r) = exitflag;
        
    end
    PDF = exp(-actions);
    PDF = PDF/trapz(Xgrid,PDF);

end