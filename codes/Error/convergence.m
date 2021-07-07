function [] = convergence(H,ERR,GERR,n,p)

switch(n)
    case(1)
        figure(100)
        hh = H;
        loglog(hh,ERR,'r','LineWidth',2)
        hold on
        grid on
        loglog(hh,hh.^(p+1),'r--','LineWidth',2)
        axis tight
        legend('||u-u_h||_{L^2(\Omega)}','h^{p+1}')
        
    case(2)
        figure(101)
        hh = H;
        loglog(hh,GERR,'b','LineWidth',2)
        hold on
        grid on
        loglog(hh,hh.^(p),'b--','LineWidth',2)
        axis tight
        legend('||u-u_h||_{H^1(\Omega)}','h^{p}')
end

