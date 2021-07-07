% Monitored Point
x = 0;
y = 3;
snap = 30;

% solution = load('TestPlaneWave_RefLev_1_snap_1.txt');
solution = load('Test_Paper_DGST_2_RefLev_1_snap_1.txt');

%select only vertex values (not all quadrature nodes) 
idSelected =  knnsearch(solution(:,1:2),[x y],'Distance','euclidean','K',1);

dis_y = zeros(snap,1);
dis_x = zeros(snap,1);

vel_y = zeros(snap,1);
for i = 1 : snap
    fname_read = ['Test_Paper_DGST_2_RefLev_1_snap_',num2str(i),'.txt'];

    disp(fname_read);
    solution = load(fname_read);  
    solutionSel = solution(idSelected,:);
    val = solutionSel(:,3:4); %dis
    dis_x(i) = val(1); 
    dis_y(i) = val(2); 
end

figure(1); hold on;
time = [0.001: 0.001 : 0.03];
plot([0 time*1000],[0; dis_y*1000]);
xlabel('time [ms]');
ylabel('u_y  [mm]')
figure(3); hold on;
time = [0.001: 0.001 : 0.03];
plot([0 time*1000],[0; dis_x*1000]);
xlabel('time [ms]');
ylabel('u_x  [mm]')

figure(2); hold on;
vel_y = [0; diff(dis_y)/0.001];
plot(time*1000,vel_y);
xlabel('time [ms]');
ylabel('v_y  [mm/s]')

