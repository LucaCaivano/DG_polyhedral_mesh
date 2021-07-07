function [fname] = write_vtk(fname, x, y, z, val1, conn, name1)



% Write vtk files for Paraview

fid = fopen(fname, 'w');
%%intestazione del file, richiesta dal formato
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Comment\n');
%%tipo file: in questo caso ascii
fprintf(fid, 'ASCII\n');
%%tipo di dati: ci sono anche structured grid e altre
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
%%punti della mesh
fprintf(fid,'POINTS %i float\n', length(x));
for i=1:length(x)
    fprintf(fid,'%g %g %g \n', x(i), y(i), z(i));
end
%if nargin < 5
%    tri = delaunay(x,y,z)-1; %%l'indice per vtk parte da 0
%end
numpoly = length(conn);
numval = 0;
for i=1:numpoly
    numval = numval + size(conn{i},1) + 1;
end

%%stampo la connettività di griglia
fprintf(fid,'CELLS %i %i\n',numpoly, numval);
for i=1:numpoly
    numpoint = size(conn{i},1);
    fprintf(fid,'%i ',[numpoint [conn{i}-1]']);
    fprintf(fid,'\n');
end
%%tipo di celle: immaigno 5 indichi il triangolo, ma è da verificare
fprintf(fid,'CELL_TYPES %i\n', numpoly);
for i=1:numpoly
    fprintf(fid,'7\n' );
end

%%%%%%
%CELL_TYPES
% 1 = vertex
% 2 = poly_vertex
% 3 = line
% 4 = poly line
% 5 = triangle
% 6 = triangle strip
% 7 = polygon
% 8 = square
% 9 = quad
%10 = tetra
%11 = voxel
%12 = hexa

%%%%%
%%dati puntuali => uno per nodo della griglia
fprintf(fid,'POINT_DATA %i\n', length(x));

%%scalars
fprintf(fid,['SCALARS  ', name1,'_x float 1\n']);
%%lookup table non so bene cosa sia, ma ci vuole
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:length(x)
    fprintf(fid,'%g\n', val1(i,1));
end
%scalars
fprintf(fid,['SCALARS  ', name1,'_y float 1\n']);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:length(x)
    fprintf(fid,'%g\n', val1(i,2));
end

%vectors
fprintf(fid,['VECTORS  ' ,name1, '  float\n']);
for i=1:length(x)
    fprintf(fid,'%g %g %g\n',val1(i,1), val1(i,2), 0);
end



fclose(fid);




