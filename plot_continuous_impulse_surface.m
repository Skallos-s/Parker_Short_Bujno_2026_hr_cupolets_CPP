% Save directory
bin_rn_direc = "save_direct/bins_1600_rN_16";

% Whether or not to plot impulses with kick
DO_KICK = false;

% Checkerboard size
Cx = 64;
Cy = 2;

% Which dimension of impulse to plot
% 2,3,4 for X,Y,Z
DIMENSION = 4;

% Directory to impulse files
full_direc = bin_rn_direc + "/bezier_impulses/";

% Suffix for files
suffix     = "_kick.txt";

if DO_KICK
    full_direc = full_direc + "kick/bin_";
else
    full_direc = full_direc + "no_kick/bin_";
    suffix = "_no" + suffix;
end

% Surface coordinates initialization
Y = [];
Z = [];

% Loop over files
for i = 0:Cx:1599

    % Load impulse file
    IMPULSE = readmatrix(full_direc + i + suffix);
    
    % Get time and x/y/z
    T = IMPULSE(:,1);
    V = IMPULSE(:,DIMENSION);
    
    % Concatenate onto coordinates
    Y = cat(2, Y, T);
    Z = cat(2, Z, V);
end

% Surface coordinates
X = meshgrid(0:Cx:1599, 1:numel(T));

% Plot surface
figure;
surf(X,Y,Z);