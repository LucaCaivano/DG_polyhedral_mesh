function [G]= SaveSolution(Dati,femregion,u_h,v_h,time)

t = time;
nln=femregion.nln;
ne=femregion.ne;

% 1D and 2D quadrature nodes and weights 
[nodes_1D, w_1D, nodes_2D, w_2D]=quadrature(Dati.nqn);

% G = [ x | y | uh(x,y) | vh(x,y) | uex(x,y) | vex(x,y) ];
G = [];

index_shift = 0;


for ie = 1 : femregion.ne % loop over elements
%     clc;
%     disp(['Computing solution for element ', num2str(ie),'/',num2str(femregion.ne)] );
    

    BBox_ie = femregion.BBox(ie,:); 
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    
    coords_elem=femregion.coords_element{ie};
    
    local_uh1 = u_h(index);
    local_uh2 = u_h(index+femregion.ndof);
    local_vh1 = v_h(index);
    local_vh2 = v_h(index+femregion.ndof);
    
    edges = [[1:femregion.nedges(ie)]' [2:femregion.nedges(ie) 1]'];
    Tria_Del = DelaunayTri(coords_elem(:,1),coords_elem(:,2), edges);
    io = Tria_Del.inOutStatus();
    Tria = Tria_Del.Triangulation(io==1,:);
            
    for iTria = 1:size(Tria,1)
        v1 = coords_elem(Tria(iTria,1),:);
        v2 = coords_elem(Tria(iTria,2),:);
        v3 = coords_elem(Tria(iTria,3),:);
        [BJ, dummy, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
        Jdet = det(BJ);
        [dphiq,Grad]= evalshape2D(femregion, ie, pphys_2D);
        for k=1:length(w_2D) % loop over quadrature nodes
            dx = abs(Jdet).*w_2D(k);
            x=pphys_2D(k,1);
            y=pphys_2D(k,2);
%             local_exact1=eval(Dati.exact_sol_1)';
%             local_exact2=eval(Dati.exact_sol_2)';
            local_aprox1=0;
            local_aprox2=0;
            local_aprox3=0;
            local_aprox4=0;
            for s=1:nln  % reconstruct the discrete solution and its gradient at the quadrature nodes
                local_aprox1 = local_aprox1 + dphiq(k,s).*local_uh1(s);
                local_aprox2 = local_aprox2 + dphiq(k,s).*local_uh2(s);
                local_aprox3 = local_aprox3 + dphiq(k,s).*local_vh1(s);
                local_aprox4 = local_aprox4 + dphiq(k,s).*local_vh2(s);
                
            end
            G = [G; x, y, local_aprox1, local_aprox2, local_aprox3, local_aprox4];
        end
    end
    
end







