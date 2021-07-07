function f = evaluate_f(neighbour,femregion,Dati,time)

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
    
    Tria_Del           = DelaunayTri(coords_elem(:,1),coords_elem(:,2), edges);
    io                 = Tria_Del.inOutStatus();
    Tria               = Tria_Del.Triangulation(io==1,:);
    
    young_ie       = Dati.Young(femregion.E_tag(ie));                               % young modulus
    poisson_ie     = Dati.Poisson(femregion.E_tag(ie));                             % poisson coefficient
    mu_ie          = young_ie/(2*(1+poisson_ie));                                   % second Lamè coefficient evaluation
    lambda_ie      = young_ie*poisson_ie/((1-2*poisson_ie)*(1+poisson_ie)) ;        % first Lamè coefficient evaluation
   
    
    for iTria = 1:size(Tria,1)
        
        v1 = coords_elem(Tria(iTria,1),:);
        v2 = coords_elem(Tria(iTria,2),:);
        v3 = coords_elem(Tria(iTria,3),:);
        
        [BJ, BJinv, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
        Jdet=det(BJ);                       % determinant
        
        % Funzioni di Base
        [dphiq,Grad] = evalshape2D(femregion, ie, pphys_2D);
        
        for k=1:length(w_2D) % loop over 2D quadrature nodes
            
            dx = w_2D(k)*Jdet;
            x = pphys_2D(k,1);
            y = pphys_2D(k,2);
            
            phi = dphiq(k,:);
            
            F1 = eval(Dati.source_1);%.*eval(Dati.source_t);
            F2 = eval(Dati.source_2);%.*eval(Dati.source_t);
            
            for i = 1:femregion.nln % loop over scalar shape functions
                 
                f1(index(i)) = f1(index(i)) + F1 * phi(i) .*dx;
                f2(index(i)) = f2(index(i)) + F2 * phi(i) .*dx;
                
                
            end
        end
    end
    
    for iedg = 1 : neighbour.nedges(ie) % loop over faces
        
        if  (neigh_ie(iedg) == -2 || neigh_ie(iedg) == -3  ||  neigh_ie(iedg) == -4) % boundary conditions
            
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
            
            for k = 1 : nqn_1D   % loop over 1D quadrature nodes
                
                ds = meshsize(iedg)*w_1D(k);
                
                Bedge = B_edge(k,:);
                
                if (neigh_ie(iedg) == -2)
                
                for i = 1 : femregion.nln % loop over scalar shape functions
                    
                    x = pphys_1D(k,1);
                    y = pphys_1D(k,2);
                    
                    neu_1 = eval(Dati.neu_1);
                    neu_2 = eval(Dati.neu_2);
                    
                    f1(index(i)) = f1(index(i)) + Bedge(i) .* neu_1 .* ds;
                    f2(index(i)) = f2(index(i)) + Bedge(i) .* neu_2 .* ds;
                    
                end
                
                
                elseif (neigh_ie(iedg) == -3)
                
                for i = 1 : femregion.nln % loop over scalar shape functions
                    
                    x = pphys_1D(k,1);
                    y = pphys_1D(k,2);
                    
                    neu_2 = eval(Dati.neu_2);
                    
                    f2(index(i)) = f2(index(i)) + Bedge(i) .* neu_2 .* ds;
                    
                end
                
                
                
                
                
                elseif (neigh_ie(iedg) == -4)
                
                for i = 1 : femregion.nln % loop over scalar shape functions
                    
                    x = pphys_1D(k,1);
                    y = pphys_1D(k,2);
                    
                    neu_1 = eval(Dati.neu_1);
                    
                    f1(index(i)) = f1(index(i)) + Bedge(i) .* neu_1 .* ds;
                    
                end
                
                
              end
           end
        end
    end
end


f = [f1; f2];