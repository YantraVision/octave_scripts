% Main Script 
% Paper:Constant Modulus Algorithm for Peak-to-Average Power Ratio (PAPR) 
% Reduction in MIMO OFDM/A
% http://ens.ewi.tudelft.nl/pubs/seyran13spl.pdf
% Seyran Khademi and Alle-Jan van der Veen 
% April 2012
%--------------------------------------------------------------------------

clear all;
close all;
clc;

%------------------ Define Parameters and Algorithms-----------------
% the default setting is for WiMAX system with OFDM block of size 1024 and 
% 60 clusters each containing 14 subcarriers. 

Mt=2;           % Number of transmit antenna
% scaling factor to change the # of subcarriers per OFDM block
% keep the N as a factor of 2 for the sake of fft calculations
% e.g. scalfac=[1,0.5,0.25] for N=[1024,512,256] respectively.
scalfac=1;      
N=1024*scalfac;         % FFT size
Nuseful=840*scalfac;    % Used subcarriers minus guards
lGuard=ceil(93*scalfac);      % Left guard subcarriers
% total number of cluster in whole OFDM symbol
Cluster_N=60*scalfac;   
M=Cluster_N;
Nb=(Nuseful/Cluster_N);    % Number of subcarriers in one cluster
Oversampling=1; % Oversampling factor
Nsym=1;         % One cluster spans Nsym OFDM symbols,It only works with 1 here
Npilot=2*Nsym;  % Number of pilots in one cluster
Nreal=10;      % Number of realizations to simulate(# of OFDM blocks)        
% Number of clusters in one disjoint block(intiger # of clusters is allowed) 
RB_in_disjoint_block=(Nuseful/M)/Nb;
% Total number of pilot subcarriers in N subcarriers
Npilot_tot=Npilot*(Cluster_N); 
% Instead of doing PAPR reduction process for whole OFDM symbols,just those
% of them who exceed the threshold can be modified 
threshold=8;%dB                
Modulation='QPSK';             % Indicates Modulation Thechnique
% Indicates Optimization Method (UnitCircle or BlockCma )
Optimization_method='BlockCma';     

%------------------------------------------------------------

% Pilot positions in total OFDM symbol
 
 if Nsym==2
 Position=ceil([5 9;1 14].*scalfac);
 else
 Position=ceil([5 9].*scalfac);   
 end
 Npilot_position=zeros(Nsym,Cluster_N*length(Position));
 
 for k=0:Cluster_N-1
 %Position of pilot symbols in a cluster of Nsym*Ncluster elements
 Npilot_position(:,(2*k+1):2*(k+1))=Position+k*(Nb); 
 end
 Npilot_position=Npilot_position+(lGuard-1);

% Generate signal constellations (64QAM, 16QAM, QPSK)
const_perm=nchoosek((1)*[1 -1 1j -1j],2);


% Pick the real and imaginary pairs and combine into one
QAM_const=const_perm(find(sign(abs(real(const_perm(:,1)))).*... %#ok
    sign(abs(imag(const_perm(:,2))))),:); %#ok
QAM_const=sum(QAM_const,2);

% Pre-allocate variables for speed
PAPR=zeros(2,Nreal);
Wrecord=zeros(Nreal,M*Mt);     
Mean=zeros(1,Nreal);
tav=0; % average processing time of the CMA algorithm
%******************************* Main Loop ********************************
for Loop=1:Nreal
    
    % Generate QPSK symbol
    % inside the () generates a random number between 1-4(for QPSK) 
    QAMs=QAM_const(ceil(length(QAM_const)*rand(Nuseful,Nsym,Mt)));
    
    % Fill upp frequency domain signal with QAMs
    SigF=zeros(N,Mt,Nsym);
    SigF(lGuard:1:(Nuseful+lGuard-1),:,:)=QAMs;
    
    % Set the pilot subcarriers to random symbols(known to receiver)
    SigF(Npilot_position,:,:)=QAM_const(ceil(length(QAM_const)...
        *rand(length(Npilot_position),Mt,Nsym)));
    
    %fill up the data matrix D
    Dprim=zeros(Nuseful,M*Mt,Nsym);
    I=zeros(M,M);
    for k=0:M-1
        I(k+1,k+1)=1;
        Q=kron(I,SigF((Nb*k)+lGuard:Nb*(k+1)+lGuard-1,:,:)); 
        Dprim=Dprim+Q;
        I(k+1,k+1)=0;
    end
    DT=zeros(N,M*Mt,Nsym);
    DT(lGuard:end-lGuard+1,:,:)=Dprim;
    
    
    % Weights (these are the spatial beamforming weights).
    % Here we assume eigen beamforming so W is a set of 
    % orthonormal vectors. Other type of beamforming could
    % be used here but eigen beamforming slightly changes the PAPR
    % property of the original OFDM signal where ZF beamformers dramatically 
    % increases the PAPR
    W=zeros(Mt,M*Mt);
    for k=0:M-1
        W(:,k*Mt+1:(k+1)*Mt)=get_orthonormal(Mt,Mt);
    end
    
    % formation of D matrix
    D=DT.';
    DNorm=norm(D,'fro');
    X=W*D; %beamformed signal
    XNorm=norm(X,'fro');
 
   % OFDM block
   F=sqrt(N).*ifft(eye(N));
   Y=X*F';
   YNorm=norm(Y,'fro');
   SigTD=Y(:);
       

   % Caculating the PAPR befor reduction
    Power1=abs(SigTD).^2;
    Pav1=mean(Power1(:));
    Pmax1=max(Power1(:));
    PAPR(1,Loop)=10*log10(Pmax1/Pav1);
    
   % formation of matrix A 
    B=D*F';
    AT=kr(B.',W); % matrix A in the paper
    ANorm=norm(AT,'fro');  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (PAPR(1,Loop)>threshold) 
% The signal goes through the PAPR reduction procedure if exceeds the 
% predefined threshold
    
      
    
    %*********** OPTIMAIZATION of PHASE SEQUENCE ************

         
     switch Optimization_method
     
     case{'none'}    
      WH=ones(M,1);   
       An=AT.';
     
     case {'BlockCma'}
          
         maxiter=50;
         alpha=Pav1;
         An=AT.';
         tic         
         [WH,Err,S]= cma_block(An,maxiter,[],alpha);           
         t=toc;
         tav=(tav+t)/2;         
         hold on
         plot(Err)  


         
    case {'UnitCircle'}     
         maxiter=50;
         alpha=Pav1;
         An=AT.';
         tic         
         [WH, Err,S] = UnitCircle(An,maxiter,[],alpha);         
         t=toc;
         tav=(tav+t)/2;
         hold on
         plot(Err)
   
     end
       %************************Optimization Finished**********************
       
       %Apply the weight coefficients to RB
       % first scale the W matrix inorder to keep the mean power constant
       WH=WH'./max(abs(W)); 
       s=WH*An;
       Transmit_Sig=s;
       
        % Calculate peak-to-average-power-ratios (PAPRs)
        
        Power2=(abs(Transmit_Sig).^2).';
        Pav2=mean(Power2(:));
        Pmax2=max(Power2(:));
        PAPR(2,Loop)=10*log10(Pmax2/Pav2);
   
       
        end
       
        fprintf('\n PAPR=%d \n',PAPR(:,Loop));
        fprintf('\n %d Realization Done\n',Loop);
end

%******************************* END of Main Loop *************************
%Plot results

figure(1)
xlabel('Iterations');
ylabel('Ojective');
title(['CMA convergence for PAPR reduction for OFDM symbol of size',num2str(N)]);
grid

figure(2)
semilogy(sort(PAPR(1,:)),1-linspace(0,1,Nreal),'-b',...
         sort(PAPR(2,:)),1-linspace(0,1,Nreal),'k')
grid
legend_text={'Original PAPR',['Precoded (',num2str(M),'RBs)']};
legend(legend_text);
xlabel('PAPR (dB)');
ylabel('CCDF');
title(['CCDF of PAPR for OFDM symbol of size',num2str(N)]);
hold on

