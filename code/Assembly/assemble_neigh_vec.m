% Authors:
% Michele Precuzzi, Luca Caivano
%--------------------------------------------------------------------

function [M]= assemble_neigh_vec(M,neight,M1,nln,n_edge)

%  INPUT:     - M: resulting matrix
%             - neight: indexes of the adjacent elements
%             - M1: matrix that collected the "neighbour component"
%             - nln: number of shape function per element
%             - n_edge: number of edges of the current element
            

for iedg=1:n_edge        % for every edge of the element
    if neight(iedg) > 0  % if the current edge is not on the boundary
        j=(neight(iedg)-1)*nln*ones(nln,1) + [1:nln]';   
        M(j) = M(j) + M1(:,iedg);
    end
end



