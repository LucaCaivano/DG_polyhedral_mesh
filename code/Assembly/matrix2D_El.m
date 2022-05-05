function [Matrices, count_2D, count_1D]= matrix2D_El(femregion,neighbour,Dati)

[nodes_1D, w_1D, nodes_2D, w_2D]=quadrature(Dati.nqn);
nqn_1D=length(w_1D);



% \int_{\Omega} (sigma(u) eps(v) dx}
V1=sparse(femregion.ndof,femregion.ndof);
V2=sparse(femregion.ndof,femregion.ndof);
V3=sparse(femregion.ndof,femregion.ndof);
V4=sparse(femregion.ndof,femregion.ndof);

% \int_{E_h} {sigma(v)} . [u]ds
IT1=sparse(femregion.ndof,femregion.ndof);
IT2=sparse(femregion.ndof,femregion.ndof);
IT3=sparse(femregion.ndof,femregion.ndof);
IT4=sparse(femregion.ndof,femregion.ndof);

% \int_{E_h} penalty  h_e^(-1) [v].[u] ds
S1=sparse(femregion.ndof,femregion.ndof);
S2=sparse(femregion.ndof,femregion.ndof);
S3=sparse(femregion.ndof,femregion.ndof);
S4=sparse(femregion.ndof,femregion.ndof);




%--------------------------------------------------------------------------------------------%

index_shift=0;
k_debug = 1;
dx_debug = [];
count_2D=0;
count_1D=0;





for ie = 1 : femregion.ne % loop over elements
    
    
    index         = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]'; %indexes of the nln degrees of freedom to work on (used instead og igloo)
    index_element = index_shift + [1:1:femregion.nedges(ie)]';  %still can't get it. Same as index but with nedges elements
    index_shift   = index_element(end);
    
    neigh_ie      = neighbour.neigh{ie};                        %indexes of the neighbour elements, counter-clockwise starting from the bottom
    neighedges_ie = neighbour.neighedges{ie};                   %indexes of the edges w.r.t. the neighbour elements
    coords_elem  = femregion.coords_element{ie};                %coordinates of the nodes of the current element
    
    [normals,meshsize] = get_normals_meshsize_faces(coords_elem);
    edges              = [[1:femregion.nedges(ie)]' [2:femregion.nedges(ie) 1]'];
    
    Tria_Del           = DelaunayTri(coords_elem(:,1),coords_elem(:,2), edges);
    
    %Tria_Del has 3 fields:
       %X: coordinates of the element's mesh points
       %Triangulation: set of indexes to indicate how it splits a square into two triangles
       %Constraints: seems to be the same as edges
    
    io                 = Tria_Del.inOutStatus();               % ???
    Tria               = Tria_Del.Triangulation(io==1,:);      % same as Tria_Del.triangulation
    
    n_qp_ie=length(w_2D)*size(Tria,1);
    count_2D = count_2D + n_qp_ie;
    
    young_ie       = Dati.Young(femregion.E_tag(ie));                               % young modulus
    poisson_ie     = Dati.Poisson(femregion.E_tag(ie));                             % poisson coefficient
    mu_ie     = young_ie/(2*(1+poisson_ie));                                        % second Lamè coefficient evaluation
    lambda_ie = young_ie*poisson_ie/((1-2*poisson_ie)*(1+poisson_ie));              % first Lamè coefficient evaluation
    
    
    
    for iTria = 1:size(Tria,1)   %loop on the subtriangles
        
        %coordinates of the subtriangle's points
        v1 = coords_elem(Tria(iTria,1),:);  
        v2 = coords_elem(Tria(iTria,2),:);
        v3 = coords_elem(Tria(iTria,3),:);
        
        [BJ, BJinv, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);  %pphys_2D contains the coordinates of the nodes of the subtriangle
        Jdet=det(BJ);                       % determinant
        
        % Basis functions 
        [dphiq,Grad] = evalshape2D(femregion, ie, pphys_2D);  %dphiq(i,j)= phi_j(pphys_2D(i))  
                                                              %Grad(i,j,k)=j-th component of the gradient of phi_k(pphys_2D(i))
        
        for k=1:length(w_2D) % loop over 2D quadrature nodes
            
            dx = w_2D(k)*Jdet;
            x = pphys_2D(k,1);
            y = pphys_2D(k,2);
            
            phi = dphiq(k,:);         % basis functions evaluated in pphys_2D(k)
            grad_x = Grad(k,1,:);     
            grad_y = Grad(k,2,:);     
                                      
            for i = 1 : femregion.nln % loop over scalar shape functions
                
                
                for j = 1 : femregion.nln % loop over scalar shape functions
                    
                    V1(index(i),index(j)) = V1(index(i),index(j)) + ((lambda_ie + 2*mu_ie)*(grad_x(j) * grad_x(i)) + mu_ie * (grad_y(j) * grad_y(i))) .*dx;
                    V2(index(i),index(j)) = V2(index(i),index(j)) + (lambda_ie * (grad_y(j) * grad_x(i)) + mu_ie * (grad_x(j) * grad_y(i))) .*dx;
                    V3(index(i),index(j)) = V3(index(i),index(j)) + (mu_ie * (grad_y(j) * grad_x(i)) + lambda_ie * (grad_x(j) * grad_y(i))) .*dx;               
                    V4(index(i),index(j)) = V4(index(i),index(j)) + ((lambda_ie + 2*mu_ie)*(grad_y(j) * grad_y(i)) + mu_ie * (grad_x(j) * grad_x(i))) .*dx;

                    
                    
                end
            end
        end
    end
    
    ITN1 = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));  % 'N' stays for 'neighbour'
    ITN2 = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));
    ITN3 = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));
    ITN4 = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));
    
    SN1 = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));
    SN4 = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));
    
    for iedg = 1:neighbour.nedges(ie)      % loop over element's edges
        
        neigedge        = neighedges_ie(iedg);   % index of the edge w.r.t. the neighbour
        ie_neigh        = neigh_ie(iedg);        % index of neighbour element
        
        n_qp_edge=length(w_1D);
        count_1D = count_1D + n_qp_edge;
        
        %%%%%%%%%%% Scaling: paper Cangiani, Georgoulis, Houston
        Cinv = femregion.area(ie)./femregion.max_kb{ie};
         
        if ie_neigh < 0 
            
             pen_coeff       = Dati.penalty_coeff.*(femregion.fem.^2).*(lambda_ie + 2 * mu_ie);  
             penalty_scaled  = pen_coeff * Cinv(iedg)*(meshsize(iedg)/femregion.area(ie));   

             young_ie_neigh        = Dati.Young(femregion.E_tag(ie));                                                          % If we are on the boundary, the element's properties are sufficient  
             poisson_ie_neigh      = Dati.Poisson(femregion.E_tag(ie));                                                        % 
             mu_ie_neigh           = young_ie_neigh/(2*(1+poisson_ie_neigh));                                                  % 
             lambda_ie_neigh       = young_ie_neigh*poisson_ie_neigh/((1-2*poisson_ie_neigh)*(1+poisson_ie_neigh));            %                      
                                                                                                                               %
        else                                                                                                                   %
                                                                                                                               %
            young_ie_neigh       = Dati.Young(femregion.E_tag(ie_neigh));                                                      %
            poisson_ie_neigh     = Dati.Poisson(femregion.E_tag(ie_neigh));                                                    %
            mu_ie_neigh          = young_ie_neigh/(2*(1+poisson_ie_neigh));                                                    %Otherwise consider also the neighbour element
            lambda_ie_neigh      = young_ie_neigh*poisson_ie_neigh/((1-2*poisson_ie_neigh)*(1+poisson_ie_neigh));       
                                                                                      
                                                                                      
            pen_coeff       = Dati.penalty_coeff.*(femregion.fem.^2).*(lambda_ie + 2 * mu_ie);
            pen_coeff_neigh = Dati.penalty_coeff.*(femregion.fem.^2).*(lambda_ie_neigh + 2 * mu_ie_neigh);
            
            Cinv_ext = femregion.area(ie_neigh) ./ femregion.max_kb{ie_neigh}(neigedge);
            
            s1 = pen_coeff * Cinv(iedg) * (meshsize(iedg)/femregion.area(ie));
            s2 = pen_coeff_neigh * Cinv_ext * (meshsize(iedg)/femregion.area(ie_neigh));
            penalty_scaled = max([s1 s2]);     
           
        end
        %%%%%%%%%%%%%%%%%%%%
        
        if iedg < neighbour.nedges(ie)   %If there are other edges of the element to explore
            p1 = coords_elem(iedg,:)'; p2 = coords_elem(iedg+1,:)';  %endpoints of the edge
            mean = 0.5*(p1+p2);                                      %middle point
            vx = p2(1)-p2(1); vy = p2(2)-p2(2); v =[vx;vy];          %direction of the edge
           % vx = p2(1)-p1(1); vy = p2(2)-p1(2); v =[vx;vy];
            v_hat = [-vy;vx];                                        %normal direction
            p3 = mean+v_hat;                                         
            v = [p1';p2';p3'];
        else                             %If this is the last edge of the element
            p1 = coords_elem(iedg,:)'; p2 = coords_elem(1,:)';       %endpoints of the edge
            mean = 0.5*(p1+p2);
           %vx = p2(1)-p2(1); vy = p2(2)-p2(2); v =[vx;vy];
            vx = p2(1)-p1(1); vy = p2(2)-p1(2); v =[vx;vy];
            v_hat = [-vy;vx];
            p3 = mean+v_hat;
            v = [p1';p2';p3'];
        end
        
        [pphys_1D] = get_jacobian_physical_points_faces(v, nodes_1D);   % quadrature nodes on the edge
        [B_edge,G_edge] = evalshape2D(femregion, ie, pphys_1D);         % shape functions and their gradients evaluated on pphys_1D
        
        if ie_neigh > 0 %ie_neigh ~= -1 && ie_neigh ~= -2
            [B_edge_neigh] = evalshape2D(femregion, neigh_ie(iedg), pphys_1D);   %If not boundary, consider also the adjacent element
        end
        
        
        for k = 1 : nqn_1D   % loop over 1D quadrature nodes
            
            ds = meshsize(iedg)*w_1D(k);     
            
            Bedge = B_edge(k,:);             %
            Gedge_x = G_edge(k,1,:);         % Evaluation on the k-th quadrature node
            Gedge_y = G_edge(k,2,:);         %
  
         
            
            for i = 1 : femregion.nln % loop over shape functions
                
                lambda_ave = 2*lambda_ie*lambda_ie_neigh/(lambda_ie + lambda_ie_neigh);
                mu_ave = 2*mu_ie*mu_ie_neigh/(mu_ie + mu_ie_neigh);
                
                aa = 0.5 * (lambda_ave + 2*mu_ave) * normals(1,iedg);     %
                ff = 0.5 * (lambda_ave + 2*mu_ave) * normals(2,iedg);     %
                bb = 0.5 * lambda_ave * normals(1,iedg);                  %
                gg = 0.5 * lambda_ave * normals(2,iedg);                  % 
                ee = 0.5 * mu_ave * normals(1,iedg);                      %
                cc = 0.5 * mu_ave * normals(2,iedg);                      %
                
                
                for j = 1 : femregion.nln % loop over scalar shape functions
                    
                    if neigh_ie(iedg) >= -1 
                        S1(index(i),index(j)) = S1(index(i),index(j)) + penalty_scaled .* Bedge(j) .* Bedge(i) .* ds;
                        S4(index(i),index(j)) = S4(index(i),index(j)) + penalty_scaled .* Bedge(j) .* Bedge(i) .* ds;
                        
                        
                    elseif neigh_ie(iedg) == -3
                        S1(index(i),index(j)) = S1(index(i),index(j)) + penalty_scaled .* Bedge(j) .* Bedge(i) .* ds;
                        
                        
                    elseif neigh_ie(iedg) == -4
                        S4(index(i),index(j)) = S4(index(i),index(j)) + penalty_scaled .* Bedge(j) .* Bedge(i) .* ds;
                        
                        
                    end
                    
                    if ie_neigh > 0 
                        
                        Bedgeneigh = B_edge_neigh(k,:);
                        
                        IT1(index(i),index(j)) = IT1(index(i),index(j)) + ( aa*Gedge_x(i).*Bedge(j) + cc*Gedge_y(i).*Bedge(j) ).* ds;
                        IT2(index(i),index(j)) = IT2(index(i),index(j)) + ( ee*Gedge_y(i).*Bedge(j) + gg*Gedge_x(i).*Bedge(j) ).* ds;
                        IT3(index(i),index(j)) = IT3(index(i),index(j)) + ( bb*Gedge_y(i).*Bedge(j) + cc*Gedge_x(i).*Bedge(j) ).* ds;
                        IT4(index(i),index(j)) = IT4(index(i),index(j)) + ( ee*Gedge_x(i).*Bedge(j) + ff*Gedge_y(i).*Bedge(j) ).* ds;
                        
                        ITN1(i,j,iedg) = ITN1(i,j,iedg) - ( aa*Gedge_x(i).*Bedgeneigh(j) + cc*Gedge_y(i).*Bedgeneigh(j) ).* ds;
                        ITN2(i,j,iedg) = ITN2(i,j,iedg) - ( ee*Gedge_y(i).*Bedgeneigh(j) + gg*Gedge_x(i).*Bedgeneigh(j) ).* ds;
                        ITN3(i,j,iedg) = ITN3(i,j,iedg) - ( bb*Gedge_y(i).*Bedgeneigh(j) + cc*Gedge_x(i).*Bedgeneigh(j) ).* ds;
                        ITN4(i,j,iedg) = ITN4(i,j,iedg) - ( ee*Gedge_x(i).*Bedgeneigh(j) + ff*Gedge_y(i).*Bedgeneigh(j) ).* ds;
                        
                        SN1(i,j,iedg) = SN1(i,j,iedg) - penalty_scaled .* Bedge(i) .* Bedgeneigh(j) .* ds;
                        SN4(i,j,iedg) = SN4(i,j,iedg) - penalty_scaled .* Bedge(i) .* Bedgeneigh(j) .* ds;
                        
                    elseif ie_neigh == -1 % boundary faces
                        
                        IT1(index(i),index(j)) = IT1(index(i),index(j)) + 2*( aa*Gedge_x(i).*Bedge(j) + cc*Gedge_y(i).*Bedge(j) ).* ds;
                        IT2(index(i),index(j)) = IT2(index(i),index(j)) + 2*( ee*Gedge_y(i).*Bedge(j) + gg*Gedge_x(i).*Bedge(j) ).* ds;
                        IT3(index(i),index(j)) = IT3(index(i),index(j)) + 2*( bb*Gedge_y(i).*Bedge(j) + cc*Gedge_x(i).*Bedge(j) ).* ds;
                        IT4(index(i),index(j)) = IT4(index(i),index(j)) + 2*( ee*Gedge_x(i).*Bedge(j) + ff*Gedge_y(i).*Bedge(j) ).* ds;
                       
                      
                    elseif ie_neigh == -3 % boundary faces
                        
                        IT1(index(i),index(j)) = IT1(index(i),index(j)) + 2*( aa*Gedge_x(i).*Bedge(j) + cc*Gedge_y(i).*Bedge(j) ).* ds;
                        IT3(index(i),index(j)) = IT3(index(i),index(j)) + 2*( bb*Gedge_y(i).*Bedge(j) + cc*Gedge_x(i).*Bedge(j) ).* ds;
                       
                        
                    elseif ie_neigh == -4 % boundary faces
                        
                        IT2(index(i),index(j)) = IT2(index(i),index(j)) + 2*( ee*Gedge_y(i).*Bedge(j) + gg*Gedge_x(i).*Bedge(j) ).* ds;
                        IT4(index(i),index(j)) = IT4(index(i),index(j)) + 2*( ee*Gedge_x(i).*Bedge(j) + ff*Gedge_y(i).*Bedge(j) ).* ds;
                   
                    end
                end
            end
        end
    end

    
    [IT1] = assemble_neigh(IT1, index, neigh_ie, ITN1, femregion.nln, neighbour.nedges(ie)); % assemble the neighbours local matrices
    [IT2] = assemble_neigh(IT2, index, neigh_ie, ITN2, femregion.nln, neighbour.nedges(ie));
    [IT3] = assemble_neigh(IT3, index, neigh_ie, ITN3, femregion.nln, neighbour.nedges(ie));
    [IT4] = assemble_neigh(IT4, index, neigh_ie, ITN4, femregion.nln, neighbour.nedges(ie));
    
    [S1] = assemble_neigh(S1, index, neigh_ie, SN1, femregion.nln, neighbour.nedges(ie));    % assemble the neighbours local
    [S4] = assemble_neigh(S4, index, neigh_ie, SN4, femregion.nln, neighbour.nedges(ie));
end






%% stiffness matrix
V = [V1, V2; V3, V4];
IT = [IT1, IT2; IT3, IT4];
S = [S1, S2; S3, S4];


Matrices = struct('V',  V , ... 
    'IT', IT, ...
    'S', S);
end

