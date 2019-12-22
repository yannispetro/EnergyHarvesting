function [Pf,shft] = get_Pf(x1,f1,parOO)
    x_failure      = parOO.x_failure;
    n_points_xmesh = parOO.n_points_xmesh;
    
    if x_failure > 0
        x_WPI = linspace(x1(1), x1(end), n_points_xmesh);
        pdf_x_WPI = interp1(x1, f1, x_WPI, 'PCHIP');

        % From here Pf
        L = 2*x_failure;

        ddx = x_WPI(2)-x_WPI(1);
        Lx = x_WPI(end)-x_WPI(1);
        xrect = linspace(-Lx/2,Lx/2,n_points_xmesh);

        rect = rectpuls(xrect,L);
        cf0 = conv(rect,pdf_x_WPI);
        [mx,ix] = max(cf0);

        shft = (n_points_xmesh - ix)*ddx;
        w = rectpuls(xrect+shft,L);

        Pf = 1 - trapz(x_WPI,pdf_x_WPI.*w);
    else
        Pf = 0;
        shft = 0;
    end
end

