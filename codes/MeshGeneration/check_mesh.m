function region = check_mesh(region)

for i = 1 %1: region.ne
    
    nedge = length(region.connectivity{i});
    
    region.coords_element{i} = region.coord(region.connectivity{i},:);
    
    max_kb{i} = zeros(nedge,1);
    
    for j = 1:nedge
       
        if j< nedge
            v1 = region.coords_element{i}(j,:); v2 = region.coords_element{i}(j+1,:);
            ch_j = j;
            ch_j_1 = j +1;
        else
            v1 = region.coords_element{i}(j,:); v2 = region.coords_element{i}(1,:);
            ch_j = j;
            ch_j_1 = 1;
        end
        for k = 1:nedge
            
            
            
            %             if k~=j && k~=(j+1)
            if k~= ch_j && k~= ch_j_1
                
%                 disp([j,k,nedge]);
                
                v3 = region.coords_element{i}(k,:);
                
                [x_tria,y_tria]=poly2cw([v1(1) v2(1) v3(1)],[v1(2) v2(2) v3(2)]);
                plot([x_tria x_tria(1)],[y_tria y_tria(1)],'b','MarkerSize',16);
                hold on;
                plot([region.coords_element{i}(end:-1:1,1); region.coords_element{i}(end,1)] ...
                     ,[region.coords_element{i}(end:-1:1,2); region.coords_element{i}(end,2)],'r','MarkerSize',16);
               
                [x1,y1] = polybool('intersection',region.coords_element{i}(end:-1:1,1),region.coords_element{i}(end:-1:1,2),x_tria,y_tria);
                
                hold on;
                plot([x1 ;x1(1)],[y1 ;y1(1)],'k','MarkerSize',16)
                
                
                area = polyarea(x_tria,y_tria);
%                 disp([polyarea(x1,y1),area]);
                
                
                if 1-any(isnan(x1)) && abs(polyarea(x1,y1)- area) < 1e-6
                    
%                     disp([max_kb{i}(j),area]);
                    close all;
                    if area>max_kb{i}(j)
                        max_kb{i}(j) = area;
                    end
                end
                
            end
        end
        
        max_kb{i}(j);
        
    end
    
end