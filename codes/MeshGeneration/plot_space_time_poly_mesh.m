function plot_space_time_poly_mesh(Dati,region, neighbour)
%figure;
% axis equal;
% axis off;
hold on;
for i = 1:size(region.coords_element,2)
    XX = region.coords_element{1,i}(:,1);
    YY = region.coords_element{1,i}(:,2);
    ZZ = 0.*region.coords_element{1,i}(:,2);
    ZZTOP = 0.*region.coords_element{1,i}(:,2) + 0.25;


    C = ones(size(XX));
    patch(XX,YY,ZZTOP, C,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.5,'Edgecolor',[0 0 0],'LineWidth',1);
end
%     plot3([XX ; XX(1)],[YY ; YY(1)],[ZZ ; ZZ(1)],'Color','k','LineWidth',1);
%     plot3([XX ; XX(1)],[YY ; YY(1)],[ZZTOP ; ZZTOP(1)],'Color','k','LineWidth',1);
%     for j = 1 : size(XX,1)
%          plot3([XX(j) XX(j)], [YY(j) YY(j)], [ZZ(j) ZZTOP(j)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-');
%     end
    %     text(mean(XX), mean(YY), num2str(i)); 
    
end
axis equal;
xlim([Dati.domain(1) Dati.domain(2)]);
ylim([Dati.domain(3) Dati.domain(4)]);

for i = 1:size(region.coords_element,2)
    X = [];
    Y = [];
    for j = 1:size(region.connectivity{i})
        % Absorbing
        if neighbour.neigh{i}(j) == -3
            if j < size(region.connectivity{i},1)
                plot([region.coords_element{i}(j,1),region.coords_element{i}(j+1,1)],[region.coords_element{i}(j,2),region.coords_element{i}(j+1,2)],'Color','g','LineWidth',1.5)
            elseif j == size(region.connectivity{i},1)
                plot([region.coords_element{i}(size(region.connectivity{i},1),1),region.coords_element{i}(1,1)],...
                    [region.coords_element{i}(size(region.connectivity{i},1),2),region.coords_element{i}(1,2)],'Color','g','LineWidth',1.5)
            end
        end
        % Neumann
        if neighbour.neigh{i}(j) == -2
            if j < size(region.connectivity{i},1)
                plot([region.coords_element{i}(j,1),region.coords_element{i}(j+1,1)],[region.coords_element{i}(j,2),region.coords_element{i}(j+1,2)],'Color','r','LineWidth',1.5)
            elseif j == size(region.connectivity{i},1)
                plot([region.coords_element{i}(size(region.connectivity{i},1),1),region.coords_element{i}(1,1)],...
                    [region.coords_element{i}(size(region.connectivity{i},1),2),region.coords_element{i}(1,2)],'Color','r','LineWidth',1.5)
            end
        end
        % Dirichlet
        if neighbour.neigh{i}(j) == -1
            if j < size(region.connectivity{i},1)
                plot([region.coords_element{i}(j,1),region.coords_element{i}(j+1,1)],[region.coords_element{i}(j,2),region.coords_element{i}(j+1,2)],'Color','b','LineWidth',1.5)
            elseif j == size(region.connectivity{i},1)
                plot([region.coords_element{i}(size(region.connectivity{i},1),1),region.coords_element{i}(1,1)],...
                    [region.coords_element{i}(size(region.connectivity{i},1),2),region.coords_element{i}(1,2)],'Color','b','LineWidth',1.5)
            end
        end
    end
end
