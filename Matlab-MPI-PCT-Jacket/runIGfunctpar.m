%function [x,fval] =  runnested(a,b,c,x0) 
%[x,fval] = fminunc(@nestedfun,x0);
%% Nested function that computes the objective function     
%    function y = nestedfun(x)
%            y = (a - b*x(1)^2 + x(1)^4/3)*x(1)^2 + x(1)*x(2) +...
%                        (-c + c*x(2)^2)*x(2)^2;     
%                            end
%                            end

% Computes the MSE for the model (eq 2)
%	A = Xb + ?
function [est,MSE] = runIGfunctpar(initVh, initVp, n, p, q, X, tau1, tau2, thread_id, options, parameter);
[est,MSE] = fminsearch(@IGfunctpar, parameter, options);
  function y = IGfunctpar(est)
    %global initVh initVp n p q X tau1 tau2
    d1=abs(est(1));
    d2=abs(est(2));
    % calculate the covariance matrices under the OU process (eq A.1)
    Vh=(d1.^tau1).*(1-d1.^(2*initVh))./(1-d1^2);
    Vp=(d2.^tau2).*(1-d2.^(2*initVp))./(1-d2^2);
    % standardize the covariance matrices to have determinant = 1 (this avoids
    % numerical issues of determinants going to infinity or zero)
    Vh=Vh./det(Vh)^(1/p);
    Vp=Vp./det(Vp)^(1/q);
    

gselect(thread_id);

% GPU computation:
gVp = gsingle(Vp);
gVh = gsingle(Vh);
gV=kron(gVp,gVh);
ginvV=inv(gV);


% compute the phylogenetic mean (eq. A.2)
gX = single(X);
gU=gones(length(X),1);
gb=(gU'*ginvV*gU)\(gU'*ginvV*gX);
% compute the phylogenetic MSE (eq. A.2)
gH=gX-gb;
gMSE=(gH'*ginvV*gH)/(n-1);

y = single(gMSE);


  end
end
