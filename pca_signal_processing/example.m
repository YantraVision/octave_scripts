%example code for explaing how pca can be used to pick up important signal sequences


N = 3400;
T = 3.4;
t = T*[0:N-1]'/N;

f1 = sin(2*pi*30*t);
f2 = 0.8*sin(2*pi*20*t);
f3 = 0.5*sin(2*pi*30*t);
f4 = 2*sin(2*pi*20*t);

f = f1 + f2 + f3 + f4; 

subplot(2,2,1)
plot(t,f1,t,f2,t,f3,t,f4)
fft_i = abs(fft(f));
subplot(2,2,2)
plot(t,fft_i);
X=[t,f1,f2,f3,f4];

%pca
[coeff, score, latent, tsquared, explained, mu] = pca(X);

%reduce the components and replot

numComp = 2;
muR = mu(1:numComp+1);
scoreR = score(:,1:numComp+1);
coeffR = coeff(1:numComp+1,1:numComp+1);
reconstR = scoreR*coeffR' + repmat(muR,3400,1);

subplot(2,2,3)
plot(t,reconstR(:,3), t, reconstR(:,2));
%plot(t,reconstR(:,2));
subplot(2,2,4)
%fft_o = abs(fft(reconstR(:,2)));
fft_o = abs(fft(reconstR(:,2) + reconstR(:,3)));
plot(t,fft_o);



