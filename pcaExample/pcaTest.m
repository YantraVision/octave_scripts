t = 0:0.1:300;
f1 = 0.2;
f2 = 0.5;
f3 = 0.3;
f4 = 150;
f5 = 300;
Fs = 1000;
wave1 = sin(2*pi*f1*t);
wave2 = sin(2*pi*f2*t);
wave3 = sin(2*pi*f3*t);
finalWave = wave1 + wave2 + wave3;
%subplot(2,2,1);
plot(t,finalWave);
title('Input Wave');
[coeff,score,latent,tsquared]=princomp(finalWave');
%printf("score value is %i\n",score);
%printf("coeff value is %i\n",coeff);
mu = mean(finalWave);
corr(score);
%printf("correlated score value is %i\n",score);
reconstructed = score * coeff' + repmat(mu, 3001,1);
%subplot(2,2,2);
%plot(t,reconstructed,"linewidth",3);
%title('Reconstructed Wave');
%printf(" value is %i\n",coeff);
