%--------------------------------------------------------------------
% PURPOSE:
%
% This routine assembles the local matrices corresponding
% to the neighbours of a given element.
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [M]= assemble_neigh(M,row,neight,M1,nln,n_edge)

%  INPUT:     - M: resulting matrix
%             - row: set of global indexes where to assembly the two matrices
%             - neight: indexes of the adjacent elements
%             - M1: matrix that collected the "neighbour component"
%             - nln: number of shape function per element
%             - n_edge: number of edges of the current element
            

for iedg=1:n_edge     % for every edge of the element
    if neight(iedg) > 0  % if the current edge is not on the boundary
        j=(neight(iedg)-1)*nln*ones(nln,1) + [1:nln]';   % e.g. [3,3,3] + [1,2,3] = [4,5,6]
        M(row,j)=M(row,j)+M1(:,:,iedg);
    end
end

