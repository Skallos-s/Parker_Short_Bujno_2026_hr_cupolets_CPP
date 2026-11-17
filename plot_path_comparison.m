% Save directory
bin_rn_direc = "save_direct/bins_1600_rN_16";
full_direc   = bin_rn_direc + "/path_comparisons/";

% Which bin to start from?
bin_start = 0:1599;

% String manipulation
if isequal(size(bin_start), [1 1])
    suffix = " " + num2str(bin_start);
else
    suffix = "s " + num2str(bin_start(1,1)) + "-" + num2str(bin_start(1,end));
end

% Create new 3D plot
figure(Name="Path Comparison Starting on Bin" + suffix);

% Loop over values of bin_start
for k = bin_start

    % Read in path data
    bezier_path = readmatrix(full_direc + "bezier_path/bin_" + k + "_bezier.txt");
    system_path = readmatrix(full_direc + "system_path/bin_" + k + "_system.txt");
    
    %Plot
    plot3(bezier_path(:,2), bezier_path(:,3), bezier_path(:,4), "w");
    hold on;
    plot3(system_path(:,2), system_path(:,3), system_path(:,4), "r");
end
hold off;