clc
clear all
close all

i_plot = 1:5;


fprintf('chain a          d          z          l          dx          Ph       Pf  \n');
for i = i_plot
%     figure(i)
    load(['computed_CHN_' num2str(i) '_xf1.4.mat'])
    
    i_out = find(computed(:,6) == 10);
    computed(i_out,:) = [];
    
    a = computed(:,1);
    d = computed(:,2);
    z = computed(:,3);
    l = computed(:,4);
    P  = -computed(:,5);
    Pf  = computed(:,6);
    dx  = computed(:,7);
    
    [Mf,If] = min(Pf);
    [M,I] = max(P(:));
    if Mf > 1e-3
        I = If;
    end
    
    ENDS(i,:) = [a(I),d(I),z(I),l(I),P(I)];
    
    fprintf('%.0f     %f & %f & %f & %f & %f & %f & %f \n', i,a(I),d(I),z(I),l(I),dx(I),P(I),Pf(I));
    
    plot1(a,d,'\alpha','\delta',[0.5,3],[0,5],1)
    plot1(a,l,'\alpha','\lambda',[0.5,3],[0,2*sqrt(5)],2)
    plot1(d,l,'\delta','\lambda',[0,5],[0,2*sqrt(5)],3)
end
figure(3)
xx = 0:0.01:5;
XX = [xx,flip(xx)];
YY = [2*sqrt(xx),2*sqrt(6)*ones(1,length(xx))];
fill(XX,YY,0.8*[1 1 1]); hold on
% plot(xx,2*sqrt(xx),'k','Linewidth',2); hold on
box off
cn = 1;
text(ENDS(cn,2)-1.3,ENDS(cn,4),['$P_h = $' num2str(ENDS(cn,5),3)],...
    'Fontsize',20,'Interpreter','latex','BackgroundColor',[1 1 1])
cn = 2;
text(ENDS(cn,2)-1.3,ENDS(cn,4),['$P_h = $' num2str(ENDS(cn,5),3)],...
    'Fontsize',20,'Interpreter','latex','BackgroundColor',[1 1 1])
cn = 3;
text(ENDS(cn,2)-1.3,ENDS(cn,4)+0.3,{'3 chains',['$P_h = $' num2str(ENDS(cn,5),3)]},...
    'Fontsize',20,'Interpreter','latex','BackgroundColor',[1 1 1])


function [] = plot1(x,y,lx,ly,lmx,lmy,i)
    figure(i)
    scatter(x,y,200,linspace(0,1,length(x)).','filled'); hold on
    set(gca,'fontsize',30)
    xlabel(lx,'fontsize',40)
    ylabel(ly,'fontsize',40)
    xlim(lmx);
    ylim(lmy);
    view(2)
end