function [neighbor]= neighbours(Dati,region,tol)

ne =region.ne;
connectivity=region.connectivity;
neigh = cell(1,ne);
neighedges = cell(1,ne);


% tagging external boundaries
% legend : -1 dirichlet, -2 neumann, -3 absorbing
left_domain  = Dati.domain(1);
right_domain = Dati.domain(2);
down_domain  = Dati.domain(3);
up_domain    = Dati.domain(4);
eps = tol*0.5;

for i=1:ne
    neigh{i}=-ones(size(connectivity{i}));
    neighedges{i}=-ones(size(connectivity{i}));
end




for i = 1:ne
    
    %   disp(['Checking element ',num2str(i),'/', num2str(ne)]);
    
    
    
    for j = 1 : length(connectivity{i})-1
        
        % Bottom boundary
        if((region.coords_element{i}(j,2) > down_domain - eps && region.coords_element{i}(j,2) < down_domain + eps) ...
                && (region.coords_element{i}(j+1,2) > down_domain - eps &&  region.coords_element{i}(j+1,2) < down_domain + eps))
            switch Dati.bc(1)
                case('D')
                    neigh{i}(j) = -1;
                    neighedges{i}(j) = -1;
                case('N')
                    neigh{i}(j) = -2;
                    neighedges{i}(j) = -2;
                case('X')
                    neigh{i}(j) = -3;
                    neighedges{i}(j) = -3;
                case('Y')
                    neigh{i}(j) = -4;
                    neighedges{i}(j) = -4;
                otherwise
                    disp('Unknown bc'); stop;
            end
        end
        
        
        
        %right
        if((region.coords_element{i}(j,1) > right_domain - eps && region.coords_element{i}(j,1) < right_domain + eps) ...
                && (region.coords_element{i}(j+1,1) > right_domain - eps &&  region.coords_element{i}(j+1,1) < right_domain + eps))
            switch Dati.bc(2)
               case('D')
                    neigh{i}(j) = -1;
                    neighedges{i}(j) = -1;
                case('N')
                    neigh{i}(j) = -2;
                    neighedges{i}(j) = -2;
                case('X')
                    neigh{i}(j) = -3;
                    neighedges{i}(j) = -3;
                case('Y')
                    neigh{i}(j) = -4;
                    neighedges{i}(j) = -4;
                otherwise
                    disp('Unknown bc'); stop;
            end
        end
        
        
        
        %top boundary
        if ((region.coords_element{i}(j,2) > up_domain - eps && region.coords_element{i}(j,2) < up_domain + eps) ...
                && (region.coords_element{i}(j+1,2) > up_domain - eps &&  region.coords_element{i}(j+1,2) < up_domain + eps))
            switch Dati.bc(3)
                case('D')
                    neigh{i}(j) = -1;
                    neighedges{i}(j) = -1;
                case('N')
                    neigh{i}(j) = -2;
                    neighedges{i}(j) = -2;
                case('X')
                    neigh{i}(j) = -3;
                    neighedges{i}(j) = -3;
                case('Y')
                    neigh{i}(j) = -4;
                    neighedges{i}(j) = -4;
                case('P')
                    x1=region.coords_element{i}(j,1);
                    x2=region.coords_element{i}(j+1,1);
                    x_mean=(x1+x2)/2;
                    if (x_mean<=1)
                       neigh{i}(j) = -4;
                       neighedges{i}(j) = -4;
                    else
                       neigh{i}(j) = -2;
                       neighedges{i}(j) = -2; 
                    end
                otherwise
                    disp('Unknown bc'); stop;
            end
        end
        
        % left boundary
        if ((region.coords_element{i}(j,1) > left_domain - eps && region.coords_element{i}(j,1) < left_domain + eps) ...
                && (region.coords_element{i}(j+1,1) > left_domain - eps &&  region.coords_element{i}(j+1,1) < left_domain + eps))
            switch Dati.bc(4)
                case('D')
                    neigh{i}(j) = -1;
                    neighedges{i}(j) = -1;
                case('N')
                    neigh{i}(j) = -2;
                    neighedges{i}(j) = -2;
                case('X')
                    neigh{i}(j) = -3;
                    neighedges{i}(j) = -3;
                case('Y')
                    neigh{i}(j) = -4;
                    neighedges{i}(j) = -4;
                otherwise
                    disp('Unknown bc'); stop;
            end
        end
        
    end
    
    
    j = length(connectivity{i});
    %down boundary
    if ((region.coords_element{i}(j,2) > down_domain - eps && region.coords_element{i}(j,2) < down_domain + eps) ...
            && (region.coords_element{i}(1,2) > down_domain - eps &&  region.coords_element{i}(1,2) < down_domain + eps))
        switch Dati.bc(1)
            case('D')
                    neigh{i}(j) = -1;
                    neighedges{i}(j) = -1;
                case('N')
                    neigh{i}(j) = -2;
                    neighedges{i}(j) = -2;
                case('X')
                    neigh{i}(j) = -3;
                    neighedges{i}(j) = -3;
                case('Y')
                    neigh{i}(j) = -4;
                    neighedges{i}(j) = -4;
            otherwise
                disp('Unknown bc'); stop;
        end
    end
    
    %right boundary
    if((region.coords_element{i}(j,1) > right_domain - eps && region.coords_element{i}(j,1) < right_domain + eps) ...
            && (region.coords_element{i}(1,1) > right_domain - eps &&  region.coords_element{i}(1,1) < right_domain + eps))
        switch Dati.bc(2)
            case('D')
                    neigh{i}(j) = -1;
                    neighedges{i}(j) = -1;
                case('N')
                    neigh{i}(j) = -2;
                    neighedges{i}(j) = -2;
                case('X')
                    neigh{i}(j) = -3;
                    neighedges{i}(j) = -3;
                case('Y')
                    neigh{i}(j) = -4;
                    neighedges{i}(j) = -4;
            otherwise
                disp('Unknown bc'); stop;
        end
    end
    
    %up boundary
    if ((region.coords_element{i}(j,2) > up_domain - eps && region.coords_element{i}(j,2) < up_domain + eps) ...
            && (region.coords_element{i}(1,2) > up_domain - eps &&  region.coords_element{i}(1,2) < up_domain + eps))
        switch Dati.bc(3)
            case('D')
                    neigh{i}(j) = -1;
                    neighedges{i}(j) = -1;
                case('N')
                    neigh{i}(j) = -2;
                    neighedges{i}(j) = -2;
                case('X')
                    neigh{i}(j) = -3;
                    neighedges{i}(j) = -3;
                case('Y')
                    neigh{i}(j) = -4;
                    neighedges{i}(j) = -4;
                case('P')
                    x1=region.coords_element{i}(j,1);
                    x2=region.coords_element{i}(1,1);
                    x_mean=(x1+x2)/2;
                    if (x_mean<=1)
                       neigh{i}(j) = -4;
                       neighedges{i}(j) = -4;
                    else
                       neigh{i}(j) = -2;
                       neighedges{i}(j) = -2; 
                    end
            otherwise
                disp('Unknown bc'); stop;
        end
    end
    
    % Left boundary
    if((region.coords_element{i}(j,1) > left_domain - eps && region.coords_element{i}(j,1) < left_domain + eps) ...
            && (region.coords_element{i}(1,1) > left_domain - eps &&  region.coords_element{i}(1,1) < left_domain + eps))
        switch Dati.bc(4)
            case('D')
                    neigh{i}(j) = -1;
                    neighedges{i}(j) = -1;
                case('N')
                    neigh{i}(j) = -2;
                    neighedges{i}(j) = -2;
                case('X')
                    neigh{i}(j) = -3;
                    neighedges{i}(j) = -3;
                case('Y')
                    neigh{i}(j) = -4;
                    neighedges{i}(j) = -4;
            otherwise
                disp('Unknown bc'); stop;
        end
    end
    
    
    
    
end




for i=1:(ne-1)
    
    
    
    edges =[];
    n_edges = length(connectivity{i});
    
    for vertices = 1:n_edges
        v(vertices)=connectivity{i}(vertices);
    end
    
    for e = 1:n_edges-1
        edges(e,:)=[v(e) v(e+1)];
    end
    edges(n_edges,:) = [v(n_edges) v(1)];
    
    for j=(i+1):ne
        edgesn =[];
        n_edgesn = length(connectivity{j});
        for verticesn = 1:n_edgesn
            vn(verticesn)=connectivity{j}(verticesn);
        end
        for e = 1:n_edgesn-1
            edgesn(e,:)=[vn(e+1) vn(e)];
        end
        edgesn(n_edgesn,:) = [vn(1) vn(n_edgesn)];
        
        for s = 1:size(edges,1)
            for t = 1:size(edgesn,1)
                if (edges(s,1) == edgesn(t,2) && edges(s,2) == edgesn(t,1))
                    neigh{i}(s)=j;
                    neigh{j}(t)=i;
                    neighedges{i}(s)=t;
                    neighedges{j}(t)=s;
                end
            end
        end
        
        for e = 1:n_edgesn-1
            edgesn(e,:)=[vn(e) vn(e+1)];
        end
        edgesn(n_edgesn,:) = [vn(n_edgesn) vn(1)];
        
        for s = 1:size(edges,1)
            for t = 1:size(edgesn,1)
                if (edges(s,1) == edgesn(t,2) && edges(s,2) == edgesn(t,1))
                    neigh{i}(s)=j;
                    neigh{j}(t)=i;
                    neighedges{i}(s)=t;
                    neighedges{j}(t)=s;
                end
            end
        end
        
    end
end
%===================================================================================
% COSTRUZIONE STRUTTURA NEIGHBOUR
%===================================================================================
neighbor.nedges = region.nedges;
neighbor.neigh = neigh;
neighbor.neighedges = neighedges;