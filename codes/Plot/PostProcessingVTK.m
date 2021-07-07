% load TestPlaneWave.mat
load Test_Paper_DGST_3.mat

connectivity = output.region.connectivity;
x = output.region.coord(:,1);
y = output.region.coord(:,2);
z = 0.*output.region.coord(:,2);

% solution = load('TestPlaneWave_RefLev_1_snap_1.txt');
solution = load('Test_Paper_DGST_3_RefLev_1_snap_1.txt');

%select only vertex values (not all quadrature nodes) 
idSelected =  knnsearch(solution(:,1:2),output.region.coord(:,1:2),'Distance','euclidean','K',1);

% scatter(solution(:,1),solution(:,2),10,solution(:,3),'filled');

for i = 1 : 3
%     fname_read = ['TestPlaneWave_RefLev_1_snap_',num2str(i),'.txt'];
    fname_read = ['Test_Paper_DGST_3_RefLev_1_snap_',num2str(i),'.txt'];

    disp(fname_read);
    solution = load(fname_read);  
    solutionSel = solution(idSelected,:);
    val = solutionSel(:,3:4); %dis
%     val = solutionSel(:,5:6); %vel

    fname = ['Snapshots_DGST_3_',num2str(i),'.vtk'];
    [fname] = write_vtk(fname, x, y, z, val, connectivity, 'dis');
end