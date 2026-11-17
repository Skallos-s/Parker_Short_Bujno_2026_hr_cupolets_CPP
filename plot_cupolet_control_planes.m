A_TC = readmatrix("bins_1600_rN_16\C001_cupolet_time_series_data.txt");
CPLT_fF = readmatrix("bins_1600_rN_16\C001_cupolet_from_impulse_data.txt");

PS2A_PT = readmatrix("bins_1600_rN_16\coding_fcn\ps2a_bin_inits.txt");

IMPULSE = readmatrix("bins_1600_rN_16\C001_cupolet_impulse_series_data.txt");

PS0 = readmatrix("control_planes\ps0_vertices.txt");
PS1 = readmatrix("control_planes\ps1_vertices.txt");
PS2 = readmatrix("control_planes\ps2_vertices.txt");

prop = 0;

T_t = A_TC(:,1); T_g = CPLT_fF(1+round(prop*end):end,1);
X_t = A_TC(:,2); X_g = CPLT_fF(1+round(prop*end):end,2);
Y_t = A_TC(:,3); Y_g = CPLT_fF(1+round(prop*end):end,3);
Z_t = A_TC(:,4); Z_g = CPLT_fF(1+round(prop*end):end,4);


clf;
plot3(X_t,Y_t,Z_t,col=[0.5,0.5,0.5])
hold on
%plot3(X_g,Y_g,Z_g,col=[1,0,0])
xlabel({'X'});
ylabel({'Y'});
zlabel({'Z'});
patch([PS0(1,1),PS0(1,1),PS0(1,1),PS0(1,1),PS0(1,1)], ...
      [PS0(1,2),PS0(1,2),PS0(2,2),PS0(2,2),PS0(1,2)], ...
      [PS0(1,3),PS0(2,3),PS0(2,3),PS0(1,3),PS0(1,3)], ...
      [0.668  0.304 0.074])
patch([PS1(1,1) PS1(1,1) PS1(2,1) PS1(2,1)], ...
      [PS1(1,2) PS1(1,2) PS1(1,2) PS1(1,2)], ...
      [PS1(1,3) PS1(2,3) PS1(2,3) PS1(1,3)], ...
      [0.773  0.617 0.125])
patch([PS2(1,1),PS2(1,1),PS2(2,1),PS2(2,1),PS2(1,1)], ...
      [PS2(1,2),PS2(1,2),PS2(1,2),PS2(2,2),PS2(1,2)], ...
      [PS2(1,3),PS2(2,3),PS2(2,3),PS2(1,3),PS2(1,3)], ...
      [0.406  0.113 0.640])
patch([PS2(3,1),PS2(3,1),PS2(4,1),PS2(4,1),PS2(3,1)], ...
      [PS2(3,2),PS2(3,2),PS2(3,2),PS2(4,2),PS2(3,2)], ...
      [PS2(3,3),PS2(4,3),PS2(4,3),PS2(3,3),PS2(1,3)], ...
      [0.316  0.578 0.207])
%x = [PS1(1,1) PS1(1,1) PS1(2,1) PS1(2,1)];
%y = [PS1(1,2) PS1(1,2) PS1(1,2) PS1(1,2)];
%z = [PS1(1,3) PS1(2,3) PS1(2,3) PS1(1,3)];



hold off

ylim([-8 2])
zlim([2.9 3.43])
%zlim([3.3 3.43])

%plot(T_t,X_t);
%hold on
%plot(T_g,X_g);
%hold off
%pause

%xlim([1.55 1.7])

%plot(IMPULSE(135838/2:end,1))
%hold on
%plot(IMPULSE(135838/2:end,1),'.')
%axis([2715 2750 -0.05 inf]);
%hold off