close all;
coeff = floor([173.38    510.67   1104.71   1755.13   2048.00   1755.13   1104.71    510.67    173.38]);

%IQ = [ 1 3 5 7 9 7 5 3 1]+1;

max_threshold = 20000;
max_val = 30155.02508450554;

scaling_val = (max_val - max_threshold);

% Cancellation signal generation
pc_signal =  (coeff*scaling_val.*IQ_set(1:9))/(2^11*max_val);

floor(pc_signal)

stem(abs(floor(IQ_set(1:9)-pc_signal)));