% Save directory
bin_rn_direc = "save_direct/bins_1600_rN_16";

% Whether or not to plot impulses with kick
DO_KICK = true;

% Color mode
% 0 for default
% 1 for gradient
% 2 for checkerboard
COLOR_MODE = 1;

% Checkerboard size
Cx = 1;
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

% Color gradient
colmap = hsv(1600);

% Create new 3D plot
figure;
plot3(0,0,0);
hold on
% Loop over files
for i = 0:Cx:1599
    if mod(i, 16) == 0
        disp(i);
    end

    % Load impulse file
    IMPULSE = readmatrix(full_direc + i + suffix);
    
    % Get time and x/y/z
    T = IMPULSE(:,1);
    X = IMPULSE(:,DIMENSION);
    
    % Plot
    switch COLOR_MODE
        case 0 % Default color
            plot3(T, X*0+i, X);

        case 1 % Rainbow gradient
            plot3(T, X*0+i, X, 'Color', colmap(i+1,:));

        case 2 % Checkerboard
            for j = 1:floor(numel(T) / Cy)
                % Checkerboard parity
                parity = mod(j + floor(i / Cx), 2);
                
                % Plot color
                col = [1,1,1] * parity;
                
                % Plot range
                R = (Cy*(j-1)+1):(Cy*j+1);
                
                plot3(T(R), X(R)*0+i, X(R), 'Color', col)
            end
    end
end
hold off