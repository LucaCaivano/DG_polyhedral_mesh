function f = evaluate_boundary_conditions(f,neighbour,femregion,Dati,time)

% INPUT:     f is the forcing term taking into account the external force
%            f1 and f2 take into account the Dirichlet condition (f1 for ux, f2 for uy)
%
% OUTPUT:    f is the total forcing term: input f + [f1;f2] 

 

[nodes_1D, w_1D, nodes_2D, w_2D]=quadrature(Dati.nqn);
nqn_1D=length(w_1D);

t = time;

f1 = sparse(femregion.ndof,1);
f2 = sparse(femregion.ndof,1);
index_shift = 0;

for ie =  1 : femregion.ne % loop over elements
    
    index         = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    index_element = index_shift + [1:1:femregion.nedges(ie)]';
    index_shift   = index_element(end);
    
    neigh_ie      = neighbour.neigh{ie};
    neighedges_ie = neighbour.neighedges{ie};
    coords_elem  = femregion.coords_element{ie};
    
    [normals,meshsize] = get_normals_meshsize_faces(coords_elem);
    edges              = [[1:femregion.nedges(ie)]' [2:femregion.nedges(ie) 1]'];
    
    
    young_ie       = Dati.Young(femregion.E_tag(ie));                               % young modulus
    poisson_ie     = Dati.Poisson(femregion.E_tag(ie));                             % poisson coefficient
    mu_ie          = young_ie/(2*(1+poisson_ie));                                   % second Lamè coefficient evaluation
    lambda_ie      = young_ie*poisson_ie/((1-2*poisson_ie)*(1+poisson_ie)) ;        % first Lamè coefficient evaluation
    
    
    
    
    for iedg = 1:neighbour.nedges(ie) % loop over faces
        
        
        %%%%%%%%%%% Scaling: paper Cangiani, Georgoulis, Houston
        Cinv = femregion.area(ie)./femregion.max_kb{ie};
        
        if (neigh_ie(iedg) == -1 || neigh_ie(iedg) == -3  ||  neigh_ie(iedg) == -4)
            
            pen_coeff       = Dati.penalty_coeff.*(femregion.fem.^2).*(lambda_ie + 2 * mu_ie);
            penalty_scaled  = pen_coeff * Cinv(iedg) * (meshsize(iedg)/femregion.area(ie));
        end
        
        %%%%%%%%%%%%%%%%%%%%
        
        if iedg < neighbour.nedges(ie)
            p1 = coords_elem(iedg,:)'; p2 = coords_elem(iedg+1,:)';
            mean = 0.5*(p1+p2);
            vx = p2(1)-p2(1); vy = p2(2)-p2(2); v =[vx;vy];
            v_hat = [-vy;vx];
            p3 = mean+v_hat;
            v = [p1';p2';p3'];
        else
            p1 = coords_elem(iedg,:)'; p2 = coords_elem(1,:)';
            mean = 0.5*(p1+p2);
            vx = p2(1)-p2(1); vy = p2(2)-p2(2); v =[vx;vy];
            v_hat = [-vy;vx];
            p3 = mean+v_hat;
            v = [p1';p2';p3'];
        end
        
        [pphys_1D] = get_jacobian_physical_points_faces(v, nodes_1D);
        [B_edge,G_edge] = evalshape2D(femregion, ie, pphys_1D);
        
        if  neigh_ie(iedg) == -1 % Dirichlet on both directions
            
            for k=1:nqn_1D   % loop over 1D quadrature nodes
                
                ds = meshsize(iedg)*w_1D(k);
                
                Bedge = B_edge(k,:);
                Gedge_x = G_edge(k,1,:);
                Gedge_y = G_edge(k,2,:);
                
                for i=1:femregion.nln % loop over scalar shape functions
                    
                    aa = 0.5 * (lambda_ie + 2*mu_ie) * normals(1,iedg);
                    ff = 0.5 * (lambda_ie + 2*mu_ie) * normals(2,iedg);
                    bb = 0.5 * lambda_ie * normals(1,iedg);
                    gg = 0.5 * lambda_ie * normals(2,iedg);
                    ee = 0.5 * mu_ie * normals(1,iedg);
                    cc = 0.5 * mu_ie * normals(2,iedg);
                    
                    
                    x=pphys_1D(k,1);
                    y=pphys_1D(k,2);
                    
                    gd1 = eval(Dati.dir_1);
                    gd2 = eval(Dati.dir_2);
                    
                    f1(index(i)) = f1(index(i)) + penalty_scaled .* Bedge(i) .* gd1 .* ds - 2*Dati.theta*( (aa*Gedge_x(i) + cc*Gedge_y(i)).*gd1 + (gg*Gedge_x(i) + ee*Gedge_y(i)).*gd2 ) .* ds;
                    f2(index(i)) = f2(index(i)) + penalty_scaled .* Bedge(i) .* gd2 .* ds - 2*Dati.theta*( (bb*Gedge_y(i) + cc*Gedge_x(i)).*gd1 + (ff*Gedge_y(i) + ee*Gedge_x(i)).*gd2 ) .* ds;
                    
                end
            end
            
            
            
        elseif  neigh_ie(iedg) == -3 % Dirichlet on x direction
            
            for k=1:nqn_1D   % loop over 1D quadrature nodes
                
                ds = meshsize(iedg)*w_1D(k);
                
                Bedge = B_edge(k,:);
                Gedge_x = G_edge(k,1,:);
                Gedge_y = G_edge(k,2,:);
                
                for i=1:femregion.nln % loop over scalar shape functions
                    
                    aa = 0.5 * (lambda_ie + 2*mu_ie) * normals(1,iedg);
                    ff = 0.5 * (lambda_ie + 2*mu_ie) * normals(2,iedg);
                    bb = 0.5 * lambda_ie * normals(1,iedg);
                    gg = 0.5 * lambda_ie * normals(2,iedg);
                    ee = 0.5 * mu_ie * normals(1,iedg);
                    cc = 0.5 * mu_ie * normals(2,iedg);
                    
                    
                    x=pphys_1D(k,1);
                    y=pphys_1D(k,2);
                    
                    gd1 = eval(Dati.dir_1);
                    %f2(index(i)) = f2(index(i)) + penalty_scaled .* Bedge(i) .* gd2 .* ds;
                     f1(index(i)) = f1(index(i)) + penalty_scaled .* Bedge(i) .* gd1 .* ds - 2*Dati.theta*( (aa*Gedge_x(i) + cc*Gedge_y(i)).*gd1 ) .* ds;
                     f2(index(i)) = f2(index(i)) - 2*Dati.theta*((bb*Gedge_y(i) + cc*Gedge_x(i)).*gd1).* ds;
                end
            end
            
            
            
            
            
        elseif  neigh_ie(iedg) == -4 % Dirichlet on y direction
            
            for k=1:nqn_1D   % loop over 1D quadrature nodes
                
                ds = meshsize(iedg)*w_1D(k);
                
                Bedge = B_edge(k,:);
                Gedge_x = G_edge(k,1,:);
                Gedge_y = G_edge(k,2,:);
                
                for i=1:femregion.nln % loop over scalar shape functions
                    
                    aa = 0.5 * (lambda_ie + 2*mu_ie) * normals(1,iedg);
                    ff = 0.5 * (lambda_ie + 2*mu_ie) * normals(2,iedg);
                    bb = 0.5 * lambda_ie * normals(1,iedg);
                    gg = 0.5 * lambda_ie * normals(2,iedg);
                    ee = 0.5 * mu_ie * normals(1,iedg);
                    cc = 0.5 * mu_ie * normals(2,iedg);
                    
                    
                    x=pphys_1D(k,1);
                    y=pphys_1D(k,2);
                    
                    gd2 = eval(Dati.dir_2);
                    
                    f1(index(i)) = f1(index(i)) -...
                        2*Dati.theta*( (gg*Gedge_x(i) + ee*Gedge_y(i)).*gd2 ) .* ds;
                    f2(index(i)) = f2(index(i)) + penalty_scaled .* Bedge(i) .* gd2 .* ds...
                        - 2*Dati.theta*((ff*Gedge_y(i) + ee*Gedge_x(i)).*gd2 ) .* ds;
               end
            end
            
            
            
            
        end
    end
end


f = f + [f1; f2];