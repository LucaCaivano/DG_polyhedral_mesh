close all

h_min=ans.hmin;
h_max=ans.hmax;
N=ans.N;
E_L2=ans.E_L2;
E_H1=ans.E_H1;
G_L2=ans.G_L2;
p=1;
n_ref=4;
figure()

% loglog((1./N), E_L2,'-o', (1./N), (1./N).^(p+1), '--', 'LineWidth',3, 'MarkerSize',10)
% title('EL2')
% legend('EL2', 'h^2')
% 
% 
% figure()
% loglog((1./N), E_H1,'-o', (1./N), (1./N).^p,'--', 'LineWidth',3,'MarkerSize',10)
% title('EH1')
% legend('EH1', 'h')
% 
% figure()
% loglog((1./N), G_L2,'-o', (1./N), (1./N).^p,'--', 'LineWidth',3,'MarkerSize',10)
% title('EH1')
% legend('G_L2', 'h')


loglog(h_min, E_L2,'-o', h_min, (h_min).^(p+1), '--', 'LineWidth',3, 'MarkerSize',10)
title('EL2')
legend('EL2', 'h^2')


figure()
loglog(h_min, E_H1,'-o', h_min, (h_min).^p,'--', 'LineWidth',3,'MarkerSize',10)
title('EH1')
legend('EH1', 'h')


 figure()
 loglog(h_min, G_L2,'-o', h_min, h_min.^p, '--', 'LineWidth',3,'MarkerSize',10)
 title('GL2')
 legend('GL2', 'h^3')

for i=1:(n_ref-1)
%    rate_of_convergence_L2(i)=-2*(log2(E_L2(i)/E_L2(i+1))/log2(N(i)/N(i+1)));
%    rate_of_convergence_H1(i)=-2*(log2(E_H1(i)/E_H1(i+1))/log2(N(i)/N(i+1)));
   rate_L2(i)=log2(E_L2(i)/E_L2(i+1));
   rate_H1(i)=log2(E_H1(i)/E_H1(i+1));
   rate_G_L2(i)=log2(G_L2(i)/G_L2(i+1));
end

