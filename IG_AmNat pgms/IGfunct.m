% Computes the MSE for the model (eq 2)%	A = Xb + ?function MSE=IGfunct(parameter)global initVh initVp n p q X tau1 tau2 invV% compute the phylogenetic mean (eq. A.2)U=ones(length(X),1);%%%size_U=size(U)%%%size_UT=size(U')%%%size_invV=size(invV)%%%size1=size(U'*invV*U)%%%size2=size(U'*invV*X)b=(U'*invV*U)\(U'*invV*X);% compute the phylogenetic MSE (eq. A.2)H=X-b;MSE=(H'*invV*H)/(n-1);