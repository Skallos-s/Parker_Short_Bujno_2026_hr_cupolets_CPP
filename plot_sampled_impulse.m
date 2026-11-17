IMPULSE = readmatrix("bins_1600_rN_16\cupolet_impulse.txt");
OFFSETS = readmatrix("bins_1600_rN_16\cupolet_offsets.txt");

S = size(IMPULSE);

% Start and end points of impulses
Starts = [];
Ends = [];

% Number of periods in cupolet time series data
PeriodCount = 50;

% Time step of cupolet
dt = 1/128;

% Impulse Dimension
% Which of x y or z is to be considered?
ImpDim = 3;

% Which impulse type to plot?
PlotImp = 7;

% Find the points
for index = 1:(S(1,1)-2)
    % Check three consecutive points
    % Impulse may incidentally be 0 during impulse
    A = IMPULSE(index+0,ImpDim);
    B = IMPULSE(index+1,ImpDim);
    C = IMPULSE(index+2,ImpDim);

    if (A == 0) && (B == 0) && (C ~= 0)
        Starts(end+1) = index+2;
    end

    if (A ~= 0) && (B == 0) && (C == 0)
        Ends(end+1) = index;
    end
end

% Don't forget ending impulse
Ends(end+1) = index+2;

% Number of impulse types
% Total number of starts is period count times impulse types
ImpTypes = numel(Starts) / PeriodCount;

% Test point
% Used for calculating offset automatically
% Shift segment from i to i+1 ontoop of each other
i = 10;
oi = IMPULSE(Starts(PlotImp) + i, 1);
oi1 = IMPULSE(Starts(PlotImp) + i + 1, 1);
odif = oi1 - oi;

figure;
plot(0,0);
title("Dimension: " + int2str(ImpDim) + " - Impulse Type: " + int2str(PlotImp))
hold on;
for o = 1:PeriodCount
    index = ImpTypes * (o - 1) + PlotImp;
    E = Ends(index);
    S = Starts(index);
    
    % Need to normalize offsets
    Offset = OFFSETS(index) / dt;

    %if IMPULSE(S+i,ImpDim) < oi
    %    Offset = (IMPULSE(S+i+1,ImpDim) - oi1) / odif;
    %end

    %if IMPULSE(S+i,ImpDim) > oi
    %    Offset = (IMPULSE(S+i,ImpDim) - oi) / odif;
    %end

    %plot(Offset + (0:E-S), IMPULSE(S:E,ImpDim))
    plot(-Offset + (0:E-S), IMPULSE(S:E,ImpDim), ".")
    %plot3(IMPULSE(S:E,1), IMPULSE(S:E,2), IMPULSE(S:E,3), ".")
end
hold off;
