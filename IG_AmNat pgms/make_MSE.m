function MSE=make_MSE(parameter)

global initVh initVp n p q X tau1 tau2

invV=parameter;

% compute the phylogenetic mean (eq. A.2)
U=ones(length(X),1);
b=(U'*invV*U)\(U'*invV*X);

% compute the phylogenetic MSE (eq. A.2)
H=X-b;
MSE=(H'*invV*H)/(n-1);