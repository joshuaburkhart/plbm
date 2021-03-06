function invV=make_Inverse(parameter)

global initVh initVp n p q X tau1 tau2

d1=abs(parameter(1));
d2=abs(parameter(2));

% calculate the covariance matrices under the OU process (eq A.1)
Vh=(d1.^tau1).*(1-d1.^(2*initVh))./(1-d1^2);
Vp=(d2.^tau2).*(1-d2.^(2*initVp))./(1-d2^2);

% standardize the covariance matrices to have determinant = 1 (this avoids
% numerical issues of determinants going to infinity or zero)
Vh=Vh./det(Vh)^(1/p);
Vp=Vp./det(Vp)^(1/q);
V=kron(Vp,Vh);

invV=V\eye(n);