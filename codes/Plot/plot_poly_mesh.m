

function plot_poly_mesh(Dati,region, neighbour)
%figure;
% axis equal;
% axis off;
hold on;
for i = 1:size(region.coords_element,2)
    XX = region.coords_element{1,i}(:,1);
    YY = region.coords_element{1,i}(:,2);
    %for j=1:size(region.coords_element{1,i}(:,1))
    %    if XX(j)<dom(1)
    %        XX(j)=dom(1);
    %    elseif XX(j)>dom(2)
    %        XX(j)=dom(2);
    %    end
    %    if YY(j)<dom(3)
    %        YY(j)=dom(3);
    %    elseif YY(j)>dom(4)
    %        YY(j)=dom(4);
    %    end
    %end
    C = ones(size(XX));
    %         patch(XX,YY,C,'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1);
    plot([XX ; XX(1)],[YY ; YY(1)],'Color','k','LineWidth',1);
%     text(mean(XX), mean(YY), num2str(i)); 
    
end
axis equal;
xlim([Dati.domain(1) Dati.domain(2)]);
ylim([Dati.domain(3) Dati.domain(4)]);

for i = 1:size(region.coords_element,2)
    X = [];
    Y = [];
    for j = 1:size(region.connectivity{i})
        % Mixed -- one Neumann one Dirichlet
        if neighbour.neigh{i}(j) == -4
            if j < size(region.connectivity{i},1)
                plot([region.coords_element{i}(j,1),region.coords_element{i}(j+1,1)],[region.coords_element{i}(j,2),region.coords_element{i}(j+1,2)],'Color','m','LineWidth',1.5)
            elseif j == size(region.connectivity{i},1)
                plot([region.coords_element{i}(size(region.connectivity{i},1),1),region.coords_element{i}(1,1)],...
                    [region.coords_element{i}(size(region.connectivity{i},1),2),region.coords_element{i}(1,2)],'Color','m','LineWidth',1.5)
            end
        end
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
