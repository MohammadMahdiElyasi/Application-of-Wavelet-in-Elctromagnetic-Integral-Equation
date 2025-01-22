%% Clear environment
clc;
clear teta f k mat ros;
close;
%% Calculating Charge Density
%%% Initializing variables
ro0 = 1; ro1 = -1; n = 512;
h = 2 * pi / n; a = 1; d = 10;
%%% Initializing pre-computations
I = 4 * (log(2) * sin(h / 2) + 2 * (sin(h / 2) * log(sin(h / 2)) - sin(h / 2)));
B = a * log(a) - a / 2 * log(2 * a ^ 2);
A = 2 * h * B;
