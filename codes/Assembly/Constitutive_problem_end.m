function[Matrices, eps_p_new_2D, eps_p_new_1D]=Constitutive_problem_end(u,Dati,eps_p_prev_2D, eps_p_prev_1D, femregion, neighbour)


%Nodes_1D= quadrature nodes in the reference 1D element
%Nodes_2D= quadrature nodes in the reference 2D element
%w_1D= weights for the 1D integration
%w_2D= weights for the 2D integration

[nodes_1D, w_1D, nodes_2D, w_2D]=quadrature(Dati.nqn);
nqn_1D=length(w_1D);

% \int_{\Omega} (T_k0(eps(u)):eps(v) dx}  using test functions for both u and v
V1=sparse(femregion.ndof,femregion.ndof);
V2=sparse(femregion.ndof,femregion.ndof);
V3=sparse(femregion.ndof,femregion.ndof);
V4=sparse(femregion.ndof,femregion.ndof);

% \int_{E_h} ({T_k0(eps(v))}.[u])ds       using test functions for both u and v
IT1=sparse(femregion.ndof,femregion.ndof);
IT2=sparse(femregion.ndof,femregion.ndof);
IT3=sparse(femregion.ndof,femregion.ndof);
IT4=sparse(femregion.ndof,femregion.ndof);

%\int_{\Omega} (T_k(u):eps(v))         using test functions for v (T_k(u) is nonlinear in u)
V_ux=sparse(femregion.ndof,1);
V_uy=sparse(femregion.ndof,1);

%\int_{E_h} ({T_k(u)}.[v])ds          using test functions for v (T_k(u) is nonlinear in u)
ITT_ux=sparse(femregion.ndof,1);
ITT_uy=sparse(femregion.ndof,1);      % the second T stays for "transpose"

index_shift=0;
qp_count_2D=0;                              % counter for the 2D nodes visited
qp_count_1D=0;                              % counter for the 1D nodes visited
eps_p_new_2D=zeros(size(eps_p_prev_2D));    % new plastic strain in the 2D nodes
eps_p_new_1D=zeros(size(eps_p_prev_1D));    % new plastic strain in the 1D nodes

n_el=0;   %number of elastic points
n_sp=0;   %number of plastic points of type 1
n_cp=0;   %number of plastic points of type 2
elastic_points=zeros(2, size(eps_p_prev_2D,2) + size(eps_p_prev_1D,2));
plastic_points=elastic_points;
critic_points=elastic_points;

for ie = 1 : femregion.ne % loop over elements
    
    index         = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]'; %indexes of the nln degrees of freedom to work on (used instead of igloo)
    index_element = index_shift + [1:1:femregion.nedges(ie)]';
    index_shift   = index_element(end);
    
    neigh_ie      = neighbour.neigh{ie};
    neighedges_ie = neighbour.neighedges{ie};
    coords_elem  = femregion.coords_element{ie};
    
    [normals,meshsize] = get_normals_meshsize_faces(coords_elem);
    edges              = [[1:femregion.nedges(ie)]' [2:femregion.nedges(ie) 1]'];
    
    Tria_Del           = DelaunayTri(coords_elem(:,1),coords_elem(:,2), edges);
    
    %Tria_Del has 3 fields:
    %X: coordinates of the element's mesh points
    %Triangulation: set of indexes to indicate how it splits a square into two triangles
    %Constraints: seems to be the same as edges
    
    io                 = Tria_Del.inOutStatus();
    Tria               = Tria_Del.Triangulation(io==1,:);      % same as Tria_Del.triangulation
    
    young_ie       = Dati.Young(femregion.E_tag(ie));                               % young modulus
    poisson_ie     = Dati.Poisson(femregion.E_tag(ie));                             % poisson coefficient
    mu_ie     = young_ie/(2*(1+poisson_ie));                                        % second Lamè coefficient evaluation
    lambda_ie = young_ie*poisson_ie/((1-2*poisson_ie)*(1+poisson_ie)) ;             % first Lamè coefficient evaluation
    
    c0_ie = Dati.Cohesion(femregion.E_tag(ie));                                        % cohesion
    f_angle_ie = Dati.Friction_angle(femregion.E_tag(ie));                             % friction angle
    K = young_ie/(3*(1-2*poisson_ie));                                                 % bulk modulus
    eta = 3*tan(f_angle_ie)/sqrt(9 + 12*tan(f_angle_ie)^2);                            % parameter for yield criterion
    c = 3*c0_ie/sqrt(9 + 12*tan(f_angle_ie)^2);                                        % parameter for yield criterion
    
    for iTria = 1:size(Tria,1)   %loop on the subtriangles
        
        %coordinates of the subtriangle's points
        v1 = coords_elem(Tria(iTria,1),:);
        v2 = coords_elem(Tria(iTria,2),:);
        v3 = coords_elem(Tria(iTria,3),:);
        
        [BJ, BJinv, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);  %pphys_2D contains the coordinates of the nodes of the subtriangle
        Jdet=det(BJ);                                                                % determinant
        
        % Basis functions (why one of them is a constant? because of modal exmapansion instead of nodal one)
        [dphiq,Grad] = evalshape2D(femregion, ie, pphys_2D);  %dphiq(i,j)= phi_j(pphys_2D(i))
        %Grad(i,j,k)=j-th component of the gradient of phi_k(pphys_2D(i))
        
        for k=1:length(w_2D) % loop over 2D quadrature nodes
            
            qp_count_2D=qp_count_2D+1;   %quadrature points index
            dx = w_2D(k)*Jdet;
            x = pphys_2D(k,1);
            y = pphys_2D(k,2);
            
            phi = dphiq(k,:);   %basis functions evaluated in pphys_2D(k)
            grad_x = Grad(k,1,:);
            grad_y = Grad(k,2,:);
            
            %Le matrici sono così divise= [xx, xy; yx, yy] e u=[ux,uy]',
            %quindi, se gli indici globali delle funzioni di base sono i1, i2, ... , inln, e
            %dovremo prendere u_i1, u_i2, ... , u_inln per ux e
            %u_(i1+ndof), u_(i2+ndof), ... , u_(inln+ndof) per uy
            
            nln=length(index);   % number of shape functions
            eps_xx=0;
            eps_xy=0;
            eps_yy=0;
            Nh=length(u)/2;      % size of u_x (or u_y)
            
            for i=1:nln
                
                eps_xx= eps_xx+u(index(i))*grad_x(:,:,i);                                            %
                eps_xy= eps_xy+(u(index(i))*grad_y(:,:,i)+u(Nh+index(i))*grad_x(:,:,i))/2;           % Assembly of total strain on the current quadrature node
                eps_yy= eps_yy+u(Nh+index(i))*grad_y(:,:,i);                                         %
                
            end
            
            eps_tot_2D=[eps_xx, eps_xy, eps_xy, eps_yy, 0]';
            [sigma_tr, p_tr, s_tr, rho_tr, n_tr]=compute_trials(eps_tot_2D, eps_p_prev_2D(:,qp_count_2D), mu_ie, lambda_ie);
            
            DP_yield = rho_tr/sqrt(2)+eta*p_tr/3-c;                   % Drucker-Prager yield function evaluation
            DP_crit = eta*p_tr/3-K*eta^2*rho_tr/(mu_ie*sqrt(2))-c;   % Second criterion evaluation (it determines if the
            % plastic point is of type 1 or 2)
            
            %CASE 1
            if(DP_yield<=0)
                
                n_el=n_el+1;
                elastic_points(:,n_el) = [x,y]';
                
                for i = 1 : femregion.nln % loop over scalar shape functions
                    
                    eps_i0=[grad_x(i), grad_y(i)/2, grad_y(i)/2, 0, 0]';
                    eps_0i=[0, grad_x(i)/2, grad_x(i)/2, grad_y(i), 0]';
                    
                    T_k = sigma_tr;
                    
                    for j = 1 : femregion.nln % loop over scalar shape functions
                        
                        V1(index(i),index(j)) = V1(index(i),index(j)) + ((lambda_ie + 2*mu_ie)*(grad_x(j) * grad_x(i)) + mu_ie * (grad_y(j) * grad_y(i))) .*dx;
                        V2(index(i),index(j)) = V2(index(i),index(j)) + (lambda_ie * (grad_y(j) * grad_x(i)) + mu_ie * (grad_x(j) * grad_y(i))) .*dx;
                        V3(index(i),index(j)) = V3(index(i),index(j)) + (mu_ie * (grad_y(j) * grad_x(i)) + lambda_ie * (grad_x(j) * grad_y(i))) .*dx;
                        V4(index(i),index(j)) = V4(index(i),index(j)) + ((lambda_ie + 2*mu_ie)*(grad_y(j) * grad_y(i)) + mu_ie * (grad_x(j) * grad_x(i))) .*dx;
                        
                    end
                    
                    V_ux(index(i))=V_ux(index(i)) + T_k'*eps_i0.*dx;
                    V_uy(index(i))=V_uy(index(i)) + T_k'*eps_0i.*dx;
                    
                end
                
                eps_p_new_2D(:,qp_count_2D)=eps_p_prev_2D(:,qp_count_2D);
                
                % CASE 2 Sysala-Valdman
                
            elseif (DP_crit<0)
                
                n_sp=n_sp+1;
                plastic_points(:,n_sp)=[x y]';
                coeff=DP_yield/(mu_ie+K*eta^2);
                T_k=sigma_tr-coeff*(mu_ie*sqrt(2)*n_tr+K*eta*[1,0,0,1,1]');
                
                for j = 1 : femregion.nln % loop over scalar shape functions
                    
                    tempxy=sqrt(2)*mu_ie*n_tr+K*eta*[1,0,0,1,1]';
                    eps_j0=[grad_x(j), grad_y(j)/2, grad_y(j)/2, 0, 0]';                                   % eps(phi_j,0)
                    tempx1=[2*grad_x(j)/3, grad_y(j)/2, grad_y(j)/2, -grad_x(j)/3, -grad_x(j)/3]';         % I_D*eps(phi_j,0)
                    tempx2=sum(n_tr.*eps_j0)*n_tr;                                                         % (n_tr x n_tr)*eps(phi_j,0)
                    tempx3=sum(tempxy.*eps_j0)*tempxy;                                                     % (tempxy x tempxy)*eps(phi,0)
                    corrector_j0=(coeff*2*sqrt(2)*mu_ie^2/rho_tr)*(tempx1-tempx2)+tempx3/(mu_ie+K*eta^2);  % T_k0*eps(phi_j,0)
                    
                    eps_0j=[0, grad_x(j)/2, grad_x(j)/2, grad_y(j), 0]';                           % eps(0,phi_j)
                    tempy1=[-grad_y(j)/3, grad_x(j)/2, grad_x(j)/2, 2*grad_y(j)/3, -grad_y(j)/3]'; % I_D*eps(0,phi_j)
                    tempy2=sum(n_tr.*eps_0j)*n_tr;                                                 % (n_tr x n_tr)*eps(0,phi_j)
                    tempy3=sum(tempxy.*eps_0j)*tempxy;                                             % (tempxy x tempxy)*eps(phi,0)
                    corrector_0j=(coeff*2*sqrt(2)*mu_ie^2/rho_tr)*(tempy1-tempy2)+tempy3/(mu_ie+K*eta^2);  % T_k0*eps(0,phi_j)
                    
                    for i = 1 : femregion.nln % loop over scalar shape functions
                        
                        eps_i0=[grad_x(i), grad_y(i)/2, grad_y(i)/2, 0, 0]';
                        eps_0i=[0, grad_x(i)/2, grad_x(i)/2, grad_y(i), 0]';
                        
                        V1(index(i),index(j)) = V1(index(i),index(j)) + ((lambda_ie + 2*mu_ie)*(grad_x(j) * grad_x(i)) + mu_ie * (grad_y(j) * grad_y(i)) - corrector_j0'*eps_i0) .*dx;
                        V2(index(i),index(j)) = V2(index(i),index(j)) + (lambda_ie * (grad_y(j) * grad_x(i)) + mu_ie * (grad_x(j) * grad_y(i)) - corrector_0j'*eps_i0) .*dx;
                        V3(index(i),index(j)) = V3(index(i),index(j)) + (mu_ie * (grad_y(j) * grad_x(i)) + lambda_ie * (grad_x(j) * grad_y(i)) - corrector_j0'*eps_0i) .*dx;
                        V4(index(i),index(j)) = V4(index(i),index(j)) + ((lambda_ie + 2*mu_ie)*(grad_y(j) * grad_y(i)) + mu_ie * (grad_x(j) * grad_x(i)) - corrector_0j'*eps_0i) .*dx;
                        
                    end
                    
                    V_ux(index(j))=V_ux(index(j)) + T_k'*eps_j0.*dx;
                    V_uy(index(j))=V_uy(index(j)) + T_k'*eps_0j.*dx;
                    
                end
                
                new_eps = eps_p_prev_2D(:,qp_count_2D)+coeff*(n_tr/sqrt(2)+eta/3*[1,0,0,1,1]');
                eps_p_new_2D(:,qp_count_2D)=new_eps;
                
                % CASE 3 Sysala-Valdman
                
            elseif (DP_crit>=0)
                
                n_cp=n_cp+1;
                critic_points(:,n_cp)=[x y]';
                
                T_k=c/eta*[1,0,0,1,1]';
                
                for i = 1 : femregion.nln % loop over scalar shape functions
                    
                    %In this case T_k0=0 -> no need for the inner for-loop
                    
                    eps_i0=[grad_x(i), grad_y(i)/2, grad_y(i)/2, 0, 0]';
                    eps_0i=[0, grad_x(i)/2, grad_x(i)/2, grad_y(i), 0]';
                    
                    V_ux(index(i))=V_ux(index(i)) + T_k'*eps_i0.*dx;
                    V_uy(index(i))=V_uy(index(i)) + T_k'*eps_0i.*dx;
                    
                end
                
                new_eps=eps_tot_2D-c/(3*K*eta)*[1,0,0,1,1]';
                eps_p_new_2D(:,qp_count_2D)=new_eps;
                
            end
        end
    end
    
    ITN1 = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));  % 'N' stays for 'neighbour'
    ITN2 = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));
    ITN3 = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));
    ITN4 = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));
    
    ITTN_ux = zeros(femregion.nln, neighbour.nedges(ie));
    ITTN_uy = zeros(femregion.nln, neighbour.nedges(ie));
    
    for iedg = 1:neighbour.nedges(ie)      % loop over element's edges
        
        neigedge        = neighedges_ie(iedg);   % index of the edge w.r.t. the neighbour
        ie_neigh        = neigh_ie(iedg);        % index of neighbour element
        
        % Scaling: paper Cangiani, Georgoulis, Houston
        Cinv = femregion.area(ie)./femregion.max_kb{ie};
        
        if ie_neigh < 0
            
            pen_coeff             = Dati.penalty_coeff.*(femregion.fem.^2).*(lambda_ie + 2 * mu_ie);
            penalty_scaled        = pen_coeff * Cinv(iedg)*(meshsize(iedg)/femregion.area(ie));
            
            young_ie_neigh       = Dati.Young(femregion.E_tag(ie));                                                          % If we are on the boundary, the element's properties are sufficient
            poisson_ie_neigh     = Dati.Poisson(femregion.E_tag(ie));                                                        %
            mu_ie_neigh          = young_ie_neigh/(2*(1+poisson_ie_neigh));                                                  %
            lambda_ie_neigh      = young_ie_neigh*poisson_ie_neigh/((1-2*poisson_ie_neigh)*(1+poisson_ie_neigh));            %
            
        else                                                                                                                  %
            
            young_ie_neigh       = Dati.Young(femregion.E_tag(ie_neigh));                                                     %
            poisson_ie_neigh     = Dati.Poisson(femregion.E_tag(ie_neigh));                                                   %
            mu_ie_neigh          = young_ie_neigh/(2*(1+poisson_ie_neigh));                                                   % Otherwise consider also the neighbour element
            lambda_ie_neigh      = young_ie_neigh*poisson_ie_neigh/((1-2*poisson_ie_neigh)*(1+poisson_ie_neigh));
            
            pen_coeff       = Dati.penalty_coeff.*(femregion.fem.^2).*(lambda_ie + 2 * mu_ie);
            pen_coeff_neigh = Dati.penalty_coeff.*(femregion.fem.^2).*(lambda_ie_neigh + 2 * mu_ie_neigh);
            
            Cinv_ext = femregion.area(ie_neigh) ./ femregion.max_kb{ie_neigh}(neigedge);
            
            s1 = pen_coeff * Cinv(iedg) * (meshsize(iedg)/femregion.area(ie));
            s2 = pen_coeff_neigh * Cinv_ext * (meshsize(iedg)/femregion.area(ie_neigh));
            penalty_scaled = max([s1 s2]);
            
        end
        
        if iedg < neighbour.nedges(ie)   %If there are other edges of the element to explore
            p1 = coords_elem(iedg,:)'; p2 = coords_elem(iedg+1,:)';  %endpoints of the edge
            mean = 0.5*(p1+p2);                                      %middle point
            vx = p2(1)-p2(1); vy = p2(2)-p2(2); v =[vx;vy];          %direction of the edge?
            v_hat = [-vy;vx];                                        %normal direction
            p3 = mean+v_hat;
            v = [p1';p2';p3'];
        else                                                         %If this is the last edge of the element
            p1 = coords_elem(iedg,:)'; p2 = coords_elem(1,:)';       %endpoints of the edge
            mean = 0.5*(p1+p2);
            vx = p2(1)-p1(1); vy = p2(2)-p1(2); v =[vx;vy];
            v_hat = [-vy;vx];
            p3 = mean+v_hat;
            v = [p1';p2';p3'];
        end
        
        [pphys_1D] = get_jacobian_physical_points_faces(v, nodes_1D);   % quadrature nodes on the edge
        [B_edge,G_edge] = evalshape2D(femregion, ie, pphys_1D);         % shape functions and their gradients evaluated on pphys_1D
        
        if ie_neigh > 0
            [B_edge_neigh, G_edge_neigh] = evalshape2D(femregion, neigh_ie(iedg), pphys_1D);   %If not boundary, consider also the adjacent element
        end
        
        for k = 1 : nqn_1D                   % loop over 1D quadrature nodes
            
            qp_count_1D=qp_count_1D+1;
            ds = meshsize(iedg)*w_1D(k);
            
            Bedge = B_edge(k,:);             % Evaluation on the k-th quadrature node
            Gedge_x = G_edge(k,1,:);         % of the shape functions of ie and their
            Gedge_y = G_edge(k,2,:);         % gradients
            
            %Gedge_neigh_x = G_edge_neigh(k,1,:);     % Evaluation on the k-th quadrature node
            %Gedge_neigh_y = G_edge_neigh(k,2,:);     % of the shape functions gradient of the neighbour
            
            nln=length(index);   % number of shape functions
            eps_xx_edge=0;
            eps_xy_edge=0;
            eps_yy_edge=0;
            Nh=length(u)/2;      % size of u_x (or u_y)
            
            for i=1:nln                                                                                     % Assembly of total strain on 1D quadrature points
                
                eps_xx_edge= eps_xx_edge+u(index(i))*Gedge_x(:,:,i);
                eps_xy_edge= eps_xy_edge+(u(index(i))*Gedge_y(:,:,i)+u(Nh+index(i))*Gedge_x(:,:,i))/2.0;
                eps_yy_edge= eps_yy_edge+u(Nh+index(i))*Gedge_y(:,:,i);
                
            end
            
            eps_tot_1D=[eps_xx_edge, eps_xy_edge, eps_xy_edge, eps_yy_edge, 0]';
            [sigma_tr, p_tr, s_tr, rho_tr, n_tr]=compute_trials(eps_tot_1D, eps_p_prev_1D(:,qp_count_1D), mu_ie, lambda_ie);
            s_tr_mat=[s_tr(1), s_tr(2), 0; s_tr(3), s_tr(4), 0; 0, 0, s_tr(5)];
            
            DP_yield=norm(s_tr_mat, 2)/sqrt(2)+eta*p_tr/3-c;                     % yield function evaluation
            DP_crit=eta*p_tr/3-K*eta^2*rho_tr/(mu_ie*sqrt(2))-c;               % second criterion evaluation
            
            lambda_ave = 2*lambda_ie*lambda_ie_neigh/(lambda_ie + lambda_ie_neigh);
            mu_ave = 2*mu_ie*mu_ie_neigh/(mu_ie + mu_ie_neigh);
            
            aa = 0.5 * (lambda_ave + 2*mu_ave) * normals(1,iedg);     %
            ff = 0.5 * (lambda_ave + 2*mu_ave) * normals(2,iedg);     %
            bb = 0.5 * lambda_ave * normals(1,iedg);                  %
            gg = 0.5 * lambda_ave * normals(2,iedg);                  %
            ee = 0.5 * mu_ave * normals(1,iedg);                      %
            cc = 0.5 * mu_ave * normals(2,iedg);                      %
            
            % CASE 1 Sysala-Valdman
            
            if(DP_yield<=0)
                
                n_el = n_el+1;
                elastic_points(:,n_el)=[x y]';
                T_k=sigma_tr;
                
                for i = 1 : femregion.nln % loop over shape functions
                    
                    for j = 1 : femregion.nln % loop over scalar shape functions
                        
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
                    
                    if ie_neigh > 0
                        
                        ITT_ux(index(i))=ITT_ux(index(i)) + 0.5*Bedge(i)*(normals(1,iedg)*T_k(1)+normals(2,iedg)*T_k(2)).*ds;
                        ITT_uy(index(i))=ITT_uy(index(i)) + 0.5*Bedge(i)*(normals(1,iedg)*T_k(3)+normals(2,iedg)*T_k(4)).*ds;
                        
                        ITTN_ux(i, iedg)=ITTN_ux(i, iedg) - 0.5*Bedgeneigh(i)*(normals(1,iedg)*T_k(1)+normals(2,iedg)*T_k(2)).*ds;
                        ITTN_uy(i, iedg)=ITTN_uy(i, iedg) - 0.5*Bedgeneigh(i)*(normals(1,iedg)*T_k(3)+normals(2,iedg)*T_k(4)).*ds;
                        
                    elseif ie_neigh == -1
                        
                        ITT_ux(index(i))=ITT_ux(index(i)) + Bedge(i)*(normals(1,iedg)*T_k(1)+normals(2,iedg)*T_k(2)).*ds;
                        ITT_uy(index(i))=ITT_uy(index(i)) + Bedge(i)*(normals(1,iedg)*T_k(3)+normals(2,iedg)*T_k(4)).*ds;
                        
                    elseif ie_neigh == -3
                        
                        ITT_ux(index(i))=ITT_ux(index(i)) + Bedge(i)*(normals(1,iedg)*T_k(1)+normals(2,iedg)*T_k(2)).*ds;
                        
                    elseif ie_neigh == -4
                        
                        ITT_uy(index(i))=ITT_uy(index(i)) + Bedge(i)*(normals(1,iedg)*T_k(3)+normals(2,iedg)*T_k(4)).*ds;
                        
                    end
                    
                end
                
                eps_p_new_1D(:,qp_count_1D)=eps_p_prev_1D(:,qp_count_1D);
                
                % CASE 2 Sysala-Valdman
                
            elseif DP_crit<0
                
                n_sp= n_sp+1;
                plastic_points(:,n_sp)=[x y]';
                coeff=DP_yield/(mu_ie+K*eta^2);
                T_k=sigma_tr-coeff*(mu_ie*sqrt(2)*n_tr+K*eta*[1,0,0,1,1]');
                
                for i = 1 : femregion.nln % loop over scalar shape functions
                    
                    tempxy=sqrt(2)*mu_ie*n_tr+K*eta*[1,0,0,1,1]';
                    
                    eps_i0=[grad_x(i), grad_y(i)/2, grad_y(i)/2, 0, 0]';                           % eps(phi_i,0)
                    tempx1=[2*grad_x(i)/3, grad_y(i)/2, grad_y(i)/2, -grad_x(i)/3, -grad_x(i)/3]'; % I_D*eps(phi_i,0)
                    tempx2=sum(n_tr.*eps_i0)*n_tr;                                                 % (n_tr x n_tr)*eps(phi_i,0)
                    tempx3=sum(tempxy.*eps_i0)*tempxy;                                             % (tempxy x tempxy)*eps(phi_i,0)
                    corrector_i0=(coeff*2*sqrt(2)*mu_ie^2/rho_tr)*(tempx1-tempx2)+tempx3/(mu_ie+K*eta^2);  % T_k0*eps(phi_i,0)
                    
                    eps_0i=[0, grad_x(i)/2, grad_x(i)/2, grad_y(i), 0]';                           % eps(0,phi_i)
                    tempy1=[-grad_y(i)/3, grad_x(i)/2, grad_x(i)/2, 2*grad_y(i)/3, -grad_y(i)/3]'; % I_D*eps(0,phi_i)
                    tempy2=sum(n_tr.*eps_0i)*n_tr;                                                 % (n_tr x n_tr)*eps(0,phi_i)
                    tempy3=sum(tempxy.*eps_0i)*tempxy;                                             % (tempxy x tempxy)*eps(0,phi_i)
                    corrector_0i=(coeff*2*sqrt(2)*mu_ie^2/rho_tr)*(tempy1-tempy2)+tempy3/(mu_ie+K*eta^2);  % T_k0*eps(0,phi_i)
                    
                    for j = 1 : femregion.nln % loop over scalar shape functions
                        
                        phi_j0_x_n=[Bedge(j)*normals(1,iedg), Bedge(j)*normals(2,iedg), 0, 0, 0]';    % tensor product between (phi_j, 0) and n+
                        phi_0j_x_n=[0, 0, Bedge(j)*normals(1,iedg), Bedge(j)*normals(2,iedg), 0]';    % tensor product between (0, phi_j) and n+
                        
                        if ie_neigh > 0 %ie_neigh ~= -1 && ie_neigh ~= -2 (internal faces)
                            
                            Bedgeneigh = B_edge_neigh(k,:);
                            
                            phi_neigh_j0_x_n=[Bedgeneigh(j)*normals(1,iedg), Bedgeneigh(j)*normals(2,iedg), 0, 0, 0]';
                            phi_neigh_0j_x_n=[0, 0, Bedgeneigh(j)*normals(1,iedg), Bedgeneigh(j)*normals(2,iedg), 0]';
                            
                            IT1(index(i),index(j)) = IT1(index(i),index(j)) + ( aa*Gedge_x(i).*Bedge(j) + cc*Gedge_y(i).*Bedge(j) - phi_j0_x_n'*corrector_i0/2).* ds;
                            IT2(index(i),index(j)) = IT2(index(i),index(j)) + ( ee*Gedge_y(i).*Bedge(j) + gg*Gedge_x(i).*Bedge(j) - phi_0j_x_n'*corrector_i0/2).* ds;
                            IT3(index(i),index(j)) = IT3(index(i),index(j)) + ( bb*Gedge_y(i).*Bedge(j) + cc*Gedge_x(i).*Bedge(j) - phi_j0_x_n'*corrector_0i/2).* ds;
                            IT4(index(i),index(j)) = IT4(index(i),index(j)) + ( ee*Gedge_x(i).*Bedge(j) + ff*Gedge_y(i).*Bedge(j) - phi_0j_x_n'*corrector_0i/2).* ds;
                            
                            ITN1(i,j,iedg) = ITN1(i,j,iedg) - ( aa*Gedge_x(i).*Bedgeneigh(j) + cc*Gedge_y(i).*Bedgeneigh(j) - phi_neigh_j0_x_n'*corrector_i0/2).* ds;
                            ITN2(i,j,iedg) = ITN2(i,j,iedg) - ( ee*Gedge_y(i).*Bedgeneigh(j) + gg*Gedge_x(i).*Bedgeneigh(j) - phi_neigh_0j_x_n'*corrector_i0/2).* ds;
                            ITN3(i,j,iedg) = ITN3(i,j,iedg) - ( bb*Gedge_y(i).*Bedgeneigh(j) + cc*Gedge_x(i).*Bedgeneigh(j) - phi_neigh_j0_x_n'*corrector_0i/2).* ds;
                            ITN4(i,j,iedg) = ITN4(i,j,iedg) - ( ee*Gedge_x(i).*Bedgeneigh(j) + ff*Gedge_y(i).*Bedgeneigh(j) - phi_neigh_0j_x_n'*corrector_0i/2).* ds;
                            
                        elseif ie_neigh == -1 % boundary faces
                            
                            IT1(index(i),index(j)) = IT1(index(i),index(j)) + ( 2*( aa*Gedge_x(i).*Bedge(j) + cc*Gedge_y(i).*Bedge(j)) - phi_j0_x_n'*corrector_i0).* ds;
                            IT2(index(i),index(j)) = IT2(index(i),index(j)) + ( 2*( ee*Gedge_y(i).*Bedge(j) + gg*Gedge_x(i).*Bedge(i)) - phi_0j_x_n'*corrector_i0).* ds;
                            IT3(index(i),index(j)) = IT3(index(i),index(j)) + ( 2*( bb*Gedge_y(i).*Bedge(j) + cc*Gedge_x(i).*Bedge(j)) - phi_j0_x_n'*corrector_0i).* ds;
                            IT4(index(i),index(j)) = IT4(index(i),index(j)) + ( 2*( ee*Gedge_x(i).*Bedge(j) + ff*Gedge_y(i).*Bedge(j)) - phi_0j_x_n'*corrector_0i).* ds;
                            
                        elseif ie_neigh == -3 % boundary faces
                            
                            IT1(index(i),index(j)) = IT1(index(i),index(j)) + ( 2*( aa*Gedge_x(i).*Bedge(j) + cc*Gedge_y(i).*Bedge(j)) - phi_j0_x_n'*corrector_i0).* ds;
                            IT3(index(i),index(j)) = IT3(index(i),index(j)) + ( 2*( bb*Gedge_y(i).*Bedge(j) + cc*Gedge_x(i).*Bedge(j)) - phi_j0_x_n'*corrector_0i).* ds;
                            
                        elseif ie_neigh == -4 % boundary faces
                            
                            IT2(index(i),index(j)) = IT2(index(i),index(j)) + ( 2*( ee*Gedge_y(i).*Bedge(j) + gg*Gedge_x(i).*Bedge(i)) - phi_0j_x_n'*corrector_i0).* ds;
                            IT4(index(i),index(j)) = IT4(index(i),index(j)) + ( 2*( ee*Gedge_x(i).*Bedge(j) + ff*Gedge_y(i).*Bedge(j)) - phi_0j_x_n'*corrector_0i).* ds;
                            
                        end
                    end
                    
                    if ie_neigh > 0
                        
                        ITT_ux(index(i))=ITT_ux(index(i)) + 0.5*Bedge(i)*(normals(1,iedg)*T_k(1)+normals(2,iedg)*T_k(2)).*ds;
                        ITT_uy(index(i))=ITT_uy(index(i)) + 0.5*Bedge(i)*(normals(1,iedg)*T_k(3)+normals(2,iedg)*T_k(4)).*ds;
                        
                        ITTN_ux(i, iedg)=ITTN_ux(i, iedg) - 0.5*Bedgeneigh(i)*(normals(1,iedg)*T_k(1)+normals(2,iedg)*T_k(2)).*ds;
                        ITTN_uy(i, iedg)=ITTN_uy(i, iedg) - 0.5*Bedgeneigh(i)*(normals(1,iedg)*T_k(3)+normals(2,iedg)*T_k(4)).*ds;
                        
                    elseif ie_neigh == -1
                        
                        ITT_ux(index(i))=ITT_ux(index(i)) + Bedge(i)*(normals(1,iedg)*T_k(1)+normals(2,iedg)*T_k(2)).*ds;
                        ITT_uy(index(i))=ITT_uy(index(i)) + Bedge(i)*(normals(1,iedg)*T_k(3)+normals(2,iedg)*T_k(4)).*ds;
                        
                        
                    elseif ie_neigh == -3
                        
                        ITT_ux(index(i))=ITT_ux(index(i)) + Bedge(i)*(normals(1,iedg)*T_k(1)+normals(2,iedg)*T_k(2)).*ds;
                        
                        
                    elseif ie_neigh == -4
                        
                        ITT_uy(index(i))=ITT_uy(index(i)) + Bedge(i)*(normals(1,iedg)*T_k(3)+normals(2,iedg)*T_k(4)).*ds;
                        
                    end
                    
                end
                
                eps_new=eps_p_prev_1D(:,qp_count_1D)+DP_yield/(mu_ie+K*eta^2)*(n_tr/sqrt(2)+eta/3*[1,0,0,1,1]');
                eps_p_new_1D(:,qp_count_1D)=eps_new;
                
                % CASE 3 Sysala-Valdman
                
            elseif DP_crit>=0
                
                n_cp= n_cp+1;
                critic_points(:,n_cp)=[x y]';
                T_k=c/eta*[1,0,0,1,1]';
                
                for i = 1 : femregion.nln % loop over shape functions
                    
                    if ie_neigh > 0
                        
                        Bedgeneigh = B_edge_neigh(k,:);
                        
                        ITT_ux(index(i))=ITT_ux(index(i)) + 0.5*Bedge(i)*(normals(1,iedg)*T_k(1)+normals(2,iedg)*T_k(2)).*ds;
                        ITT_uy(index(i))=ITT_uy(index(i)) + 0.5*Bedge(i)*(normals(1,iedg)*T_k(3)+normals(2,iedg)*T_k(4)).*ds;
                        
                        ITTN_ux(i, iedg)=ITTN_ux(i, iedg) - 0.5*Bedgeneigh(i)*(normals(1,iedg)*T_k(1)+normals(2,iedg)*T_k(2)).*ds;
                        ITTN_uy(i, iedg)=ITTN_uy(i, iedg) - 0.5*Bedgeneigh(i)*(normals(1,iedg)*T_k(3)+normals(2,iedg)*T_k(4)).*ds;
                        
                    elseif ie_neigh == -1
                        
                        ITT_ux(index(i))=ITT_ux(index(i)) + Bedge(i)*(normals(1,iedg)*T_k(1)+normals(2,iedg)*T_k(2)).*ds;
                        ITT_uy(index(i))=ITT_uy(index(i)) + Bedge(i)*(normals(1,iedg)*T_k(3)+normals(2,iedg)*T_k(4)).*ds;
                        
                    elseif ie_neigh == -3
                        
                        ITT_ux(index(i))=ITT_ux(index(i)) + Bedge(i)*(normals(1,iedg)*T_k(1)+normals(2,iedg)*T_k(2)).*ds;
                        
                    elseif ie_neigh == -4
                        
                        ITT_uy(index(i))=ITT_uy(index(i)) + Bedge(i)*(normals(1,iedg)*T_k(3)+normals(2,iedg)*T_k(4)).*ds;
                        
                    end
                end
                
                new_eps=eps_tot_1D-c/(3*K*eta)*[1,0,0,1,1]';
                eps_p_new_1D(:,qp_count_1D)=new_eps;
                
            end
        end
    end
    
    [IT1] = assemble_neigh(IT1, index, neigh_ie, ITN1, femregion.nln, neighbour.nedges(ie)); % assemble the neighbours local matrices
    [IT2] = assemble_neigh(IT2, index, neigh_ie, ITN2, femregion.nln, neighbour.nedges(ie));
    [IT3] = assemble_neigh(IT3, index, neigh_ie, ITN3, femregion.nln, neighbour.nedges(ie));
    [IT4] = assemble_neigh(IT4, index, neigh_ie, ITN4, femregion.nln, neighbour.nedges(ie));
    
    [ITT_ux] = assemble_neigh_vec(ITT_ux, neigh_ie, ITTN_ux, femregion.nln, neighbour.nedges(ie));
    [ITT_uy] = assemble_neigh_vec(ITT_uy, neigh_ie, ITTN_uy, femregion.nln, neighbour.nedges(ie));
    
end

ITT = [IT1, IT2; IT3, IT4]';
V  = [V1, V2; V3, V4];

ITT_u = [ITT_ux; ITT_uy];
V_u  = [V_ux; V_uy];


% disp('N° elastic points')
% n_el
% disp('N° plastic (simple) points')
% n_sp
% disp('N° plastic points')
% n_cp

Matrices = struct('V',  V , ...
    'ITT', ITT, ...
    'V_u', V_u, ...
    'ITT_u', ITT_u);

%Plot of elastic and plastic points
figure()
scatter(elastic_points(1, 1:n_el), elastic_points(2, 1:n_el), 'filled', 'b')
hold on
scatter(plastic_points(1, 1:n_sp), plastic_points(2, 1:n_sp), 'filled', 'r')
scatter(critic_points(1, 1:n_cp), critic_points(2, 1:n_cp), 'filled', 'g')
legend({'Elastic points', 'Plastic points', 'Critic points'})

% fprintf('  Proportion of plastic points: %6.1e, ',n_sp/(n_el+n_sp+n_cp)*100);
% fprintf('  Proportion of critic points: %6.1e, ',n_cp/(n_el+n_sp+n_cp)*100);

end