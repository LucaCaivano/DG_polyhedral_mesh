function [output] = MAIN_QN(testcase,varargin)


% NOT TESTED!



% load dati structure
Dati = dati(testcase);
Dati.nqn = 2*Dati.fem + 1; %number of quadrature GL nodes

if (strcmp(Dati.type,'verification'))
    
    disp(['Verication test case ...',testcase]);
    N = varargin{1};
    % infos about the mesh
    disp(['Number of Polygonal Elements ... ', num2str(N)]);
    
    if size(varargin,2) == 1
        refinement_step = 4;
    else
        refinement_step = varargin{2};
    end
    % infos about the mesh
    disp(['solution computed for ', num2str(refinement_step),' refinement level']);
    
else
    
    % load external mesh
    filemate = Dati.MeshName;
    load([filemate,'.mat']);
    % infos about the mesh
    disp(['Number of Polygonal Elements ... ', num2str(N)]);
    refinement_step = 1;
end

Newton_tol=Dati.tol;
Newton_super_tol=Newton_tol/100;
Newton_max_it=Dati.max_it;

for rl = 1 : refinement_step
    
    if (strcmp(Dati.type,'verification'))
        
        % only for square domains -- cartesian mesh
        h = (Dati.domain(2)-Dati.domain(1))/N;
        [X1,Y1] = meshgrid([h/2 : h : Dati.domain(2)], [h/2 : h : Dati.domain(4)]);  % middle points of the mesh elements
        P = [[X1(:)], [Y1(:)]];
        Nel = size(P,1);
        %Square with square mesh
        [region] = generate_mesh(Dati,Dati.domain,@MbbDomain,Nel,100, P);  %If varargin{1}=2 there are 4 subsquares of [0,1]^2
        
        %Square with polytopic mesh
        %[region] = generate_mesh(Dati,Dati.domain,@MbbDomain,Nel,100);
        
        %Corona
        %[region] = generate_mesh_corona(Dati,Dati.domain,@CircleDomain,Nel,100);
        
        %Plate
        %[region] = generate_mesh(Dati,Dati.domain,@MichellDomain,Nel,100);
        
    end
    
    %region's fields:
    %nedges: nedges[i] contains the number of edges of the i-th element
    %BBox: bounding box
    %ne: number of elements
    %coord: set of points' coordinates
    %coords_element: same as coord but in groups, one for each element
    %connectivity: a set of ne vectors in which we find the enumeration
    % of the internal and boundary nodes (same order as
    % coord)
    %area: set of ne values representing the elements' area
    %max_kb: used only for penality calibration
    
    hmin = min(sqrt((region.BBox(:,1)-region.BBox(:,2)).^2 + (region.BBox(:,3)-region.BBox(:,4)).^2));
    hmax = max(sqrt((region.BBox(:,1)-region.BBox(:,2)).^2 + (region.BBox(:,3)-region.BBox(:,4)).^2));
    disp(['Minimum mesh size ...', num2str(hmin)]);
    disp(['Maximum mesh size ...', num2str(hmax)]);
    
    % connectivity
    disp('Making connectivity ... ');
    [neighbour] = neighbours(Dati,region,hmin);
    disp('Done');
    
    
    %neighbour has 3 fields:
    %nedges: same as region.nedges
    %neigh: neigh{i} contains the list of the indexes of the adjacent elements (counter clockwise starting from the bottom)
    %neighedges: neghbour's edges enumeration
    
    % femregion
    disp('Making femregion ... ');
    [femregion] = create_dof(Dati,region);
    
    %femregion has several fields:
    %fem: degree of the polynomials involved
    %nedges: same as neighbour.nedges
    %nln: number of degrees of freedom in a single element
    %ndof: number of total degrees of freedom
    %ne: number of elements
    %nqn: number of quadrature nodes per element
    %coord: same as region.coord (set of mesh points' coords)
    %BBox: same as region.BBox
    %coords_element: same as region.coords_element
    %connectivity: same as region.connectivity
    %area: same a region.area
    %max_kb: same as region.max_kb
    %E_tag: used to evaluate rho, vs and vp
    
    if (strcmp(Dati.type,'verification'))
        femregion.E_tag(1:femregion.ne) = 1;
    else
        femregion.E_tag = E_tag;
    end
    disp('Done');
    
    % plot mesh with boundary conditions
    disp('Plotting the ploygonal mesh');
    if (strcmp(Dati.type,'verification'))
        subplot(1,4,rl)
        plot_poly_mesh(Dati,region, neighbour);
        
    else
        plot_poly_mesh(Dati,region, neighbour);
    end
    disp('Done');
    
    % saving output
    output.hmin(rl) = hmin;
    output.hmax(rl) = hmax;
    output.region(rl) = region;
    output.neighbour(rl) = neighbour;
    output.femregion(rl) = femregion;
    output.Dati(rl) = Dati;
    
    
    %building of  stiffness and mass matrix
    disp('Making Mass and Stiffness matrices');
    [Matrices, count_2D, count_1D] = matrix2D_El(femregion,neighbour,Dati);
    V_el=Matrices.V;
    IT_el=Matrices.IT;
    S=Matrices.S;
    K_el=V_el - Dati.theta*IT_el - IT_el' + S;
    disp('Done');
    
    eps_p_prev_2D=zeros(5, count_2D);
    eps_p_prev_1D=zeros(5, count_1D);
    
    % for QUASI-NEWTON
    rho_0 = 4.5e8;
    H = 1/rho_0*eye(size(K_el,1));      % inverse of the initial guess
    
    dim=2*femregion.ndof;
    disp('Beginning of the time loop')
    
    % loading parameters
    d_zeta=1/1000;                         % load increment
    d_zeta_min=d_zeta/1300;                % minimum increment
    d_zeta_old=d_zeta;
    zeta=0;                                % current load factor
    zeta_old=zeta;
    zeta_max=1;                            % maximal value of the load factor (full D-condition)
    
    f = evaluate_f(neighbour,femregion,Dati,0.0);
    Dirichlet_f=evaluate_boundary_conditions(zeros(dim,1),neighbour,femregion,Dati,0);  %rhs assuming the full boundary condition
    final_f = f + Dirichlet_f;
    u=zeros(dim,1);                        % solution
    initial_f = d_zeta*final_f;
    u_it=K_el\initial_f;                   % initial guess for the first loading step (elastic solution)
    u_el = K_el\final_f;
    u_old=-u_it;                           % solution at previous loading step (like an arificial u at step -1)
    du=zeros(dim,1);                       % increment at every Newton iteration
    
    
    %Plot the elastic solution
    v2 = zeros(size(u_el));
    [G] = SaveSolution(Dati,femregion,u_el,v2,0);
    filename = [testcase,'_RefLev_', '0','_elastic_','0','.txt'];
    fid = fopen(filename,'w');
    for i = 1 : size(G,1)
        fprintf(fid,'%8.7e,%8.7e,%8.7e,%8.7e,%8.7e,%8.7e', G(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
    %filename = [testcase,'_RefLev_', num2str(rl),'_time_',num2str(t),'.txt'];
    filename = [testcase,'_RefLev_', '0','_elastic_','0','.txt']; %(by me)
    sol = load(filename);
    % plot_poly_solution(output.Dati(end),output.region(end),sol,fig_num)
    figure()
    scatter(sol(:,1),sol(:,2),20,sol(:,3),'filled');colorbar
    title('u_x^{el}')
    figure()
    scatter(sol(:,1),sol(:,2),20,sol(:,4),'filled');colorbar
    title('u_y^{el}')
    figure()
    scatter(sol(:, 1), sol(:, 2), 10, sqrt(sol(:, 3).^2+sol(:, 4).^2), 'filled'), colorbar
    title('u_{tot}^{el}')
    figure()
    scatter(sol(:, 1)+sol(:, 3), sol(:, 2)+sol(:, 4), 10, sqrt(sol(:, 3).^2+sol(:, 4).^2), 'filled'), colorbar
    title('elastic strain')
    
    
    %Plot of initial solution
    v2 = zeros(size(u_it));
    [G] = SaveSolution(Dati,femregion,u_it,v2,0);
    filename = [testcase,'_RefLev_', '0','_elastic_','0','.txt'];
    fid = fopen(filename,'w');
    for i = 1 : size(G,1)
        fprintf(fid,'%8.7e,%8.7e,%8.7e,%8.7e,%8.7e,%8.7e', G(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
    filename = [testcase,'_RefLev_', '0','_elastic_','0','.txt']; %(by me)
    sol = load(filename);
    % plot_poly_solution(output.Dati(end),output.region(end),sol,fig_num)
    figure()
    scatter(sol(:,1),sol(:,2),20,sol(:,3),'filled');colorbar
    title('u_x^{el}')
    figure()
    scatter(sol(:,1),sol(:,2),20,sol(:,4),'filled');colorbar
    title('u_y^{el}')
    figure()
    scatter(sol(:, 1), sol(:, 2), 10, sqrt(sol(:, 3).^2+sol(:, 4).^2), 'filled'), colorbar
    title('u_{tot}^{el}')
    figure()
    scatter(sol(:, 1)+sol(:, 3), sol(:, 2)+sol(:, 4), 10, sqrt(sol(:, 3).^2+sol(:, 4).^2), 'filled'), colorbar
    title('elastic strain')
    
    
    n_newton=1;
    n_converged=0;
    n_non_converged=0;
    n_super_converged=0;
    
    load_hist = zeros(1,20);
    crit_hist = zeros(1,20);
    it_hist = zeros(1,20);
    
    
    disp('Beginning loading')
    
    while (zeta<zeta_max)
        
        
        zeta=zeta_old+d_zeta;          % total load we are going to consider
        f_step = zeta*final_f;
        criterion=Newton_tol+1;
        it=0;
        
        disp('Beginning Quasi-Newton')
        while (it<Newton_max_it)
            
            if (it == 0)
                [Matrices_pl]=Constitutive_problem(u_it,Dati,eps_p_prev_2D,eps_p_prev_1D,femregion,neighbour);   %compute the nonlinear parts of the stiffness matrix K_k
                %and the tangential stiffness matrix K_k0
                
                V_pl = Matrices_pl.V;
                ITT_pl = Matrices_pl.ITT;
                V_u = Matrices_pl.V_u;
                ITT_u = Matrices_pl.ITT_u;
                
                K_k0 = V_pl - ITT_pl +S;
                K_ku = V_u - ITT_u + S*u_it;
                
                r_k=f_step-K_ku;
                
            end
            
            du = -H*r_k;            % note that H ~ (-K_k0)^-1
            u_new = u + du;
            
            [Matrices_pl]=Constitutive_problem(u_new,Dati,eps_p_prev_2D,eps_p_prev_1D,femregion,neighbour);
            V_pl = Matrices_pl.V;
            ITT_pl = Matrices_pl.ITT;
            V_u = Matrices_pl.V_u;
            ITT_u = Matrices_pl.ITT_u;
            
            K_k0 = V_pl - ITT_pl +S;
            K_ku = V_u - ITT_u + S*u_it;
            
            r_k_new = f_step-K_ku;
            y_k = r_k_new - r_k;
            
            H = (eye(size(H,1))-du*y_k'/(y_k'*du))*H*(eye(size(H,1))-y_k*du'/(y_k'*du)) + du*du'/(y_k'*du);
            
            r_k = r_k_new;
            
            % stopping criterion
            q1 = sqrt( du'*K_el*du );
            q2 = sqrt( u_it'*K_el*u_it );
            q3 = sqrt( u_new'*K_el*u_new );
            criterion = q1/(q2+q3);
            
            u_it=u_new;
            
            if (criterion<Newton_super_tol)
                break;
            end
            
            it=it+1;
            
        end
        
        load_hist(n_newton) = zeta;
        crit_hist(n_newton) = criterion;
        it_hist(n_newton) = it;
        n_newton=n_newton+1;
        
        %Plot of u_it after newton
        v2 = zeros(size(u_it));
        [G] = SaveSolution(Dati,femregion,u_it,v2,0);
        filename = [testcase,'_RefLev_', '0','_elastic_','0','.txt'];
        fid = fopen(filename,'w');
        for i = 1 : size(G,1)
            fprintf(fid,'%8.7e,%8.7e,%8.7e,%8.7e,%8.7e,%8.7e', G(i,:));
            fprintf(fid,'\n');
        end
        fclose(fid);
        %filename = [testcase,'_RefLev_', num2str(rl),'_time_',num2str(t),'.txt'];
        filename = [testcase,'_RefLev_', '0','_elastic_','0','.txt']; %(by me)
        sol = load(filename);
        %plot_poly_solution(output.Dati(end),output.region(end),sol,fig_num)
        figure()
        scatter(sol(:,1),sol(:,2),20,sol(:,3),'filled');colorbar
        title('u_x^{it}')
        figure()
        scatter(sol(:,1),sol(:,2),20,sol(:,4),'filled');colorbar
        title('u_y^{it}')
        figure()
        scatter(sol(:, 1), sol(:, 2), 10, sqrt(sol(:, 3).^2+sol(:, 4).^2), 'filled'), colorbar
        title('u_{tot}^{it}')
        figure()
        scatter(sol(:, 1)+sol(:, 3), sol(:, 2)+sol(:, 4), 10, sqrt(sol(:, 3).^2+sol(:, 4).^2), 'filled'), colorbar
        title('elastic strain it')
        figure()
        subplot(1, 2, 1)
        scatter(sol(:,1),sol(:,2),20,sol(:,3),'filled');colorbar
        title('ux')
        subplot(1, 2, 2)
        scatter(sol(:,1),sol(:,2),20,sol(:,4),'filled');colorbar
        title('uy')
        pause(2)
        
        
        if (criterion <= Newton_tol)
            n_converged=n_converged+1;
            disp('Converged!')
            u_old=u;
            u=0.5*u_it;   %for SV
            %u=u_it;
            [~, eps_p_prev_2D, eps_p_prev_1D]=Constitutive_problem_end(u,Dati,eps_p_prev_2D,eps_p_prev_1D,femregion,neighbour);      % update strains
            zeta_old=zeta;
            d_zeta_old=d_zeta;
            
            if (criterion<=Newton_super_tol && n_converged>=7)
                disp('Super converged!')
                n_super_converged=n_super_converged+1;
                d_zeta=d_zeta*2;
            end
            
        else
            disp('Not converged')
            n_non_converged=n_non_converged+1;
            d_zeta=d_zeta/2;
        end
        
        % initialization for the next iteration
        u_it=d_zeta*(u-u_old)/d_zeta_old+u;  % let u_delta= u-u_old (difference between solutions at current and previous time steps)
        % then U_it= U+alpha*U_delta
        % where alpha =  1/2 if convergenve was not achieved
        %                1   if Newton_super_tol<criterion<Newton_tol
        %                2   if criterion<=Newton_super_tol
        
        if d_zeta<d_zeta_min
            warning('Too small load increments.')
            break
        end
        
    end
    
    %double the elements
    N = N*2;
end

%Plot final solution
v2 = zeros(size(u));
[G] = SaveSolution(Dati,femregion,u,v2,0);
%             filename = [testcase,'_RefLev_', num2str(rl),'_time_',num2str(t),'.txt'];
filename = [testcase,'_RefLev_', '0','_snap_','0','.txt'];
fid = fopen(filename,'w');
for i = 1 : size(G,1)
    fprintf(fid,'%8.7e,%8.7e,%8.7e,%8.7e,%8.7e,%8.7e', G(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
filename = [testcase,'_RefLev_', '0','_snap_','0','.txt']; %(by me)
sol = load(filename);
% plot_poly_solution(output.Dati(end),output.region(end),sol,fig_num)
figure()
scatter(sol(:,1),sol(:,2),20,sol(:,3),'filled');colorbar
title('u_x^{pl}')
figure()
scatter(sol(:,1),sol(:,2),20,sol(:,4),'filled');colorbar
title('u_y^{pl}')
figure()
scatter(sol(:, 1), sol(:, 2), 10, sqrt(sol(:, 3).^2+sol(:, 4).^2), 'filled'), colorbar
title('u_{tot}^{pl}')
figure()
scatter(sol(:, 1) + sol(:, 3), sol(:, 2) + sol(:, 4), 10, sqrt(sol(:, 3).^2+sol(:, 4).^2), 'filled'), colorbar
title('plastic strain')

end
