% Save directory
bin_rn_direc = "save_direct\bins_1600_rN_16";

% Whether or not to plot impulses with kick
DO_KICK = false;

% Which dimension of impulse to plot
% 2,3,4 for X,Y,Z
DIMENSION = 2;

% Directory to impulse files
full_direc = bin_rn_direc + "\bezier_impulses\";

% Suffix for files
suffix     = "_kick.txt";

if DO_KICK
    full_direc = full_direc + "kick\bin_";
else
    full_direc = full_direc + "no_kick\bin_";
    suffix = "_no" + suffix;
end

% Create new 3D plot
figure;
plot3(0,0,0);
hold on
% Loop over files
for i = 0:1599
    
    % Load impulse file
    IMPULSE = readmatrix(full_direc + i + suffix);
    
    % Get time and x/y/z
    T = IMPULSE(:,1);
    X = IMPULSE(:,DIMENSION);
    
    % Plot
    plot3(T, X*0+i, X);
end
hold off