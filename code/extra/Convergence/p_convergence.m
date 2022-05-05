p=[1, 2, 3, 4];
% res.E_H1(4)=ans.E_H1;
% res.E_L2(4)=ans.E_L2;
% res.G_L2(4)=ans.G_L2;
E_L2=res.E_L2;
E_H1=res.E_H1;
G_L2=res.G_L2;



figure()
semilogy(p, E_L2,'-o', p, exp(-p), '--', 'LineWidth',3, 'MarkerSize',10)
title('E_L2')
legend('E_L2', 'e^p')
figure()
semilogy(p, E_H1,'-o', p, exp(-p), '--', 'LineWidth',3, 'MarkerSize',10)
title('E_H1')
legend('E_H1', 'e^p')
figure()
semilogy(p, G_L2,'-o', p, exp(-p), '--', 'LineWidth',3, 'MarkerSize',10)
title('G_L2')
legend('E_L2', 'e^p')
