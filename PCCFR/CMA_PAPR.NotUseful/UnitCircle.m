
% This function implements the block version of the modified constant modulus 
% algorithm for the PAPR objective function with unimodular beamformer
% weights
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

function [W, Err, S] = UnitCircle(X,maxiter,W_init,alpha)

[M,N] = size(X);	% M = number of design parameters , N = number of subcarriers
mu =0.02 ;		% small step size (dept. on scaling of X)
mup=mu/(alpha^2);

if isempty(W_init)
  
   W_init = 2*round(rand(M,1))-1;		
end

One = 1*alpha.*ones(1,N);
W = W_init;
Err = zeros(1,maxiter);

for ii = 1:maxiter,  
	w = W;
	y = w'*X;	% should have constant modulus 1 once w is correct
	e = (conj(y).*y - One)  ; 
	ye = y .* e;
	w = w - mup * X* ye'  ;
    wnorm=abs(w);
    % inorder to get the unimodular beamformer the solution is divided by
    % its norm at each iteration
	W = w./wnorm;
	Err(ii) = norm(e);
end

S = W'*X;
