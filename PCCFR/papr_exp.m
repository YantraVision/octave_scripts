clc;
close all;
clear all;
% Signal generation and peak generation
fid  = fopen ("pdsch.out", "r", "ieee-le");
IQ = fread (fid,1024, "int16");
fclose(fid);

IQ = reshape(IQ,2,length(IQ)/2);
IQ = IQ(1,:) + i*IQ(2,:);

%Interpolation 
IQ = interp1 ( IQ, [1:0.25:length(IQ)],"spline");
IQ = 32767*IQ/max(abs(IQ));

% Impulse response
%coeff = [0.075031 0.403263 0.631495 0.861850 1.000000 0.861850 0.631495 0.403263 0.075031];
#coeff = ones(1,9);
#coeff = fir1 (8, 0.0001);
#coeff = coeff/max(coeff);
coeff = window ('gausswin', 9)';

% Caculating the PAPR befor reduction
Power1=abs(IQ).^2;
Pav1=mean(Power1(:));
Pmax1=max(Power1(:));
PAPR =10*log10(Pmax1/Pav1)

% Reference peak generation
% Peak detection
max_threshold = 30000; #use threshold
[max_val,idx] = max(IQ);

max_val = abs(IQ(451));
idx = 451;

scaling_val = (max_val - max_threshold)/max_val;

% Cancellation signal generation
pc_signal =  coeff*scaling_val.* IQ(idx-4:idx+4);


% Final compensation
IQ_out = IQ;
IQ_out(idx-4:idx+4) = IQ_out(idx-4:idx+4) - pc_signal;

subplot(2,1,1);
stem(real((IQ(idx-4:idx+4))));
hold on
stem(real((IQ_out(idx-4:idx+4))),'r*');

subplot(2,1,2);
stem(imag((IQ(idx-4:idx+4))));
hold on
stem(imag((IQ_out(idx-4:idx+4))),'r*');

figure(2)
stem(abs((IQ(idx-4:idx+4))));
hold on
stem(abs((IQ_out(idx-4:idx+4))),'r*');

figure(3)
stem(abs((IQ)));
hold on
stem(abs((IQ_out)),'r*');
plot(max_threshold*ones(1,length(IQ)),'k');

Power1=abs(IQ_out).^2;
Pav1=mean(Power1(:));
Pmax1=max(Power1(:));
PAPR =10*log10(Pmax1/Pav1)