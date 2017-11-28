close all
x = linspace(0.01,0.5,256);
X = abs(x);
x_phase = arg(x);

Y = Amp_Saleh_model(X);
hold on
plot(abs(Y),'k')
plot(X,'g')

pause 

Y_inv = X./Y;

fid  = fopen('coeff.txt','w');
fprintf(fid,'%d, %d, \n',real(Y_inv),imag(Y_inv));
fclose(fid);


plot((e.^x_phase).* abs(X.*Y.*Y_inv)/2047);
hold
plot(real(X.*Y.*Y_inv)/2047,'r+');

close all
Y_pre_distort = X.*Y_inv; % can be accesd by LUT
Y_DPD = Amp_Saleh_model(Y_pre_distort);
plot(real(Y_DPD),'r')
hold on
plot(abs(Y),'k')
plot(X,'g')
