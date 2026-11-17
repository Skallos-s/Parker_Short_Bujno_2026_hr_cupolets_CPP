IMPULSE = readmatrix("bins_1600_rN_16\C001_cupolet_impulse_series_data.txt");

S = size(IMPULSE);

% Start and end points of impulses
Starts = [];
Ends = [];

% Number of impulse types in cupolet
ImpTypes = 3;

% Impulse Dimension
% Which of x y or z is to be considered?
ImpDim = 1;

% Which impulse type to plot?
PlotImp = 3;

% Find the points
for index = 1:(S(1,1)-2)
    % Check three consecutive points
    % Impulse may incidentally be 0 during impulse
    A = IMPULSE(index+0,PlotImp);
    B = IMPULSE(index+1,PlotImp);
    C = IMPULSE(index+2,PlotImp);

    if (A == 0) && (B == 0) && (C ~= 0)
        Starts(end+1) = index+2;
    end

    if (A ~= 0) && (B == 0) && (C == 0)
        Ends(end+1) = index;
    end
end

% Don't forget ending impulse
Ends(end+1) = index+2;

% Test point
% Used for calculating offset automatically
% Shift segment from i to i+1 ontoop of each other
i = 10;
oi = IMPULSE(Starts(PlotImp) + i, 1);
oi1 = IMPULSE(Starts(PlotImp) + i + 1, 1);
odif = oi1 - oi;

figure;
hold on;
for o = 1:numel(Starts) / ImpTypes
    index = ImpTypes * (o - 1) + PlotImp;
    E = Ends(index);
    S = Starts(index);

    Offset = 0;
    if IMPULSE(S+i,ImpDim) < oi
        Offset = (IMPULSE(S+i+1,ImpDim) - oi1) / odif;
    end

    if IMPULSE(S+i,ImpDim) > oi
        Offset = (IMPULSE(S+i,ImpDim) - oi) / odif;
    end

    plot(Offset + (0:E-S), IMPULSE(S:E,ImpDim))
    plot(Offset + (0:E-S), IMPULSE(S:E,ImpDim),".")
end
hold off;
