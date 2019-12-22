clc
clear all
close all

%           a,  d,  z,   l
parOO.lb = [0.5 0.0 0.05 0];
parOO.ub = [3.0 5.0 0.20 2*sqrt(5)];
nvars = length(parOO.lb);

parOO.x_failure = 2.8/2;
parOO.Pf_acceptable = 0.001;
parOO.n_points_xmesh = 2^14+1;

% Parameters
% parS.zeta   = 0.1;
% parS.delta  = 0.2;
% parS.lambda = 2*sqrt(parS.delta);
parS.kappa  = 0.65;
% parS.gam    = 0.8;
parS.S0     = 0.05;
parS.ord    = 2;
parS.ndof   = 2;
parS.Calpha = 1;

parW.ti     = 0;
parW.tf     = 10;
parW.points = 31;
parW.Fs     = 50;
parW.h      = 13;

[g0,galpha,g1,g2] = Legendre(parW,parS.Calpha);
parW.g0 = g0;
parW.galpha = galpha;
parW.g1 = g1;
parW.g2 = g2;

parO.TolsIP  = [1e-6 1e-7 1e-6]; % OptTol StepTol ConstrTol
parO.MaxIter = 30000;
parO.MaxFun  = 90000;

n_chains = 5;

INITPs = choose_initP('Init_grid_data.mat',parOO,20,n_chains);

options = optimoptions('patternsearch','Display','iter', ...
    'FunctionTolerance',1e-5, ...
    'MeshTolerance',1e-4, ...
    'StepTolerance',eps, ...
    'InitialMeshSize',0.75, ...
    'MeshContractionFactor',0.5);

global computed
for kk = 1:n_chains
    computed = INITPs(kk,:);
    x0_PatS  = (INITPs(kk,1:4) - parOO.lb)./(parOO.ub - parOO.lb);

    obj = @(x) objPS_Pf(x, parOO, parS, parW, parO, kk);
    
    % ---- Algorithm 1
    lb0 = zeros(1,nvars);
    ub0 = ones(1,nvars);
    [best_PatS_x,PatSval,exitflag,output] = patternsearch( ...
            obj,x0_PatS,[],[],[],[],lb0,ub0,[],options);
    disp(['paternsearch - ' num2str(output.iterations) ' objective evaluations'])
    disp([best_PatS_x PatSval])% -nonlcon(best_PatS_x)])
end