function plot_poly_solution(Dati,region, solution,fig_num)
hold on;



for i = 1 : size(region.coords_element,2)
   
    
    XX = region.coords_element{1,i}(:,1);
    YY = region.coords_element{1,i}(:,2);
    
    %K = dsearchn(X,XI) returns the indices K of the closest points in X for
    %each point in XI.
    K = dsearchn(solution(:,1:2),[XX,YY]); 
    
    Cx = solution(K,3);
    Cy = solution(K,4);
    figure(fig_num);
    subplot(1,2,1)
    patch(XX,YY,Cx); hold on;
    plot([XX ; XX(1)],[YY ; YY(1)],'Color','k','LineWidth',1);
    title('x-component'); colorbar;
    subplot(1,2,2)
    patch(XX,YY,Cy); hold on;
    plot([XX ; XX(1)],[YY ; YY(1)],'Color','k','LineWidth',1);
    title('y-component'); colorbar;
    
        
end


figure(fig_num);
subplot(1,2,1)
axis equal;
xlim([Dati.domain(1) Dati.domain(2)]);
ylim([Dati.domain(3) Dati.domain(4)]);
subplot(1,2,2)
axis equal;
xlim([Dati.domain(1) Dati.domain(2)]);
ylim([Dati.domain(3) Dati.domain(4)]);
