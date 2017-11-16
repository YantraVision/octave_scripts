
% This function implements the block version of the standard constant modulus 
% algorithm CMA(2,2) (Treichler 1983) for the PAPR objective function
% Paper:Constant Modulus Algorithm for Peak-to-Average Power Ratio (PAPR) 
% Reduction in MIMO OFDM/A
% http://ens.ewi.tudelft.nl/pubs/seyran13spl.pdf
% Seyran Khademi and Alle-Jan van der Veen 
% April 2012

% W: desired estimated weight vector matrix (M * d);
% Err: improvement in objective function 
% S: recovered constant-modulus signal (as many as cols of W_init)
% W: estimated weight vector
% W_init: initial guess for the weight vector matrix (M * d);
% X: Data matrix, alpha: average power of the signal
%--------------------------------------------------------------------------

function [W,Err,S] = cma_block(X,maxiter,W_init,alpha)

Err=zeros(1,maxiter);
[M,N] = size(X);	% M = number of design parameters , N = number of subcarriers

mu =0.05 ;		% small step size (dept. on scaling of X)
mup=mu/(alpha^2); %scaling the stepsize according to the average signal power

if isempty(W_init)
   W_init = (1/sqrt(2))*(ones(M,1)+1j*ones(M,1));		% initial solution
end

One = 1*alpha.*ones(1,N);
W = W_init;


for ii = 1:maxiter,

	w = W;
	y = w'*X;	% should have constant modulus 1 once w is correct
	e = (conj(y).*y - One)  ;
	ye = y .* e;
	w = w - mup * X* ye'  ;
    W = w;
	Err(ii) = norm(e);
end

 S = W'*X;
