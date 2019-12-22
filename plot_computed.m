clc
clear all
close all

load('computed_CHN_3_xf1.4.mat')
size(computed)

computed = computed(1:end,:);

id_in = find(computed(:,6) ~= 10);
computed = computed(id_in,:);

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
format long

fprintf('a          d          z          l          dx          Ph       Pf  \n');
fprintf('%f & %f & %f & %f & %f & %f & %f \n', a(I),d(I),z(I),l(I),dx(I),P(I),Pf(I));

plot1(a,d,'\alpha','\delta',[0.5,3],[0,5])
% plot1(a,z,'\alpha','\zeta')
plot1(a,l,'\alpha','\lambda',[0.5,3],[0,2*sqrt(5)])
% plot1(d,z,'\delta','\zeta')
plot1(d,l,'\delta','\lambda',[0,5],[0,2*sqrt(5)])
% plot1(z,l,'\zeta','\lambda')

function [] = plot1(x,y,lx,ly,lmx,lmy)
    figure()
    scatter(x,y,200,[1:length(x)].','filled'); hold on
    set(gca,'fontsize',30)
    xlabel(lx,'fontsize',40)
    ylabel(ly,'fontsize',40)
    if strcmp(lx,'\delta') && strcmp(ly,'\lambda')
        xx = lmx(1):0.01:lmx(2);
        plot(xx,2*sqrt(xx))
    end
    xlim(lmx);
    ylim(lmy);
    view(2)
end