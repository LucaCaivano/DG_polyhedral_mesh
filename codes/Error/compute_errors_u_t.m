%--------------------------------------------------------------------
% PURPOSE:
%
% This routine computes the errors in different norms, i.e.,
%
% ||u-u_h||_{0,\Omega}
% ||u-u_h||_{1,\Omega}
% ||u-u_h||_{\infinity}
% 
% Author:
% Paola Antonietti
%--------------------------------------------------------------------


function [E_L2]= compute_errors(Dati,femregion,u_t,time)

nln=femregion.nln;
ne=femregion.ne;

% initialization
E_L2 = 0;
t = time;
% 1D and 2D quadrature nodes and weights 
[nodes_1D, w_1D, nodes_2D, w_2D]=quadrature(Dati.nqn);

index_shift = 0;
for ie=1:femregion.ne % loop over elements
    
    BBox_ie = femregion.BBox(ie,:); 
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    
    coords_elem=femregion.coords_element{ie};
    
    local_uh1 = u_t(index);
    local_uh2 = u_t(index+femregion.ndof);
        
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
            
            local_exact1 = eval(Dati.exact_grad1_t)';
            local_exact2 = eval(Dati.exact_grad2_t)';
            
            local_aprox1=0;
            local_aprox2=0;
            
            for s=1:nln  % reconstruct the discrete solution and its gradient at the quadrature nodes
                local_aprox1 = local_aprox1 + dphiq(k,s).*local_uh1(s);
                local_aprox2 = local_aprox2 + dphiq(k,s).*local_uh2(s);
            end

            E_L2 = E_L2 + ((local_aprox1 - local_exact1).^2).*dx + ((local_aprox2 - local_exact2).^2).*dx;

        end
    end
end

E_L2 = sqrt(E_L2);