function val = get_Ph(x2,f2,parOO,gam)
    n_points_xmesh = parOO.n_points_xmesh;

    y_WPI = linspace(x2(1), x2(end), n_points_xmesh);
    pdf_y_WPI = interp1(x2, f2, y_WPI, 'PCHIP');
    var_y_WPI = trapz(y_WPI, (y_WPI.^2).*pdf_y_WPI);
    P = gam*var_y_WPI; % Non-dimesional harvested power

    val = -P;

%         if plt == true
%             close all
%             figure(2)
%             scatter(x2,f2); hold on
%             plot(y_WPI,pdf_y_WPI); hold on
%             scatter(x1,f1); hold on
%             plot(x_WPI,pdf_x_WPI); hold on
%             plot(x_WPI,w.*max(f1),'k'); hold on
%             xlim([x1(1) x1(end)])
%             ylim([0,1.5])
%             drawnow
%         end
end

