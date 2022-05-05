function [u] = evaluate_solution(Dati,femregion,neighbour,time,ux,uy)

[nodes_1D, w_1D, nodes_2D, w_2D]=quadrature(Dati.nqn);

t=time;

u1=sparse(femregion.ndof,1);               % \int_{\Omega} f1 . v dx + boundary conditions
u2=sparse(femregion.ndof,1);               % \int_{\Omega} f2 . v dx + boundary conditions

index_shift=0;


for ie = 1 : femregion.ne % loop over elements
%     clc;
%     disp(['Computing solution for element ', num2str(ie),'/',num2str(femregion.ne)] );

    
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    index_element = index_shift + [1:1:femregion.nedges(ie)]';
    index_shift = index_element(end);
    
    neigh_ie = neighbour.neigh{ie};
    neighedges_ie = neighbour.neighedges{ie};
    coords_elem = femregion.coords_element{ie};
    
    [normals,meshsize]=get_normals_meshsize_faces(coords_elem);
    edges = [[1:femregion.nedges(ie)]' [2:femregion.nedges(ie) 1]'];
    Tria_Del = DelaunayTri(coords_elem(:,1),coords_elem(:,2), edges);
    io = Tria_Del.inOutStatus();
    Tria = Tria_Del.Triangulation(io==1,:);
    for iTria = 1:size(Tria,1)

        v1 = coords_elem(Tria(iTria,1),:);
        v2 = coords_elem(Tria(iTria,2),:);
        v3 = coords_elem(Tria(iTria,3),:);
 
        [BJ, BJinv, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
        Jdet=det(BJ);                       % determinant

        % Funzioni di Base
        [dphiq,Grad] = evalshape2D(femregion, ie, pphys_2D);
                        
        for k=1:length(w_2D) % loop over 2D quadrature nodes

            dx=w_2D(k)*Jdet;
            x=pphys_2D(k,1);
            y=pphys_2D(k,2);
            
            phi = dphiq(k,:);
            dt = Dati.dt;  
            
            U1 = eval(ux);
            U2 = eval(uy);
            
            for i=1:femregion.nln % loop over scalar shape functions
                
                u1(index(i)) = u1(index(i)) + U1 * phi(i) .*dx;
                u2(index(i)) = u2(index(i)) + U2 * phi(i) .*dx;
                
            end
            
            
        end
    end
    
    
    
    
end
u = [u1; u2];
