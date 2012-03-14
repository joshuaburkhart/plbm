% Computes the MSE for the model (eq 2)%	A = Xb + ?function MSE=IGfunct(parameter)global initVh initVp n p q X tau1 tau2d1=abs(parameter(1));d2=abs(parameter(2));% calculate the covariance matrices under the OU process (eq A.1)Vh=(d1.^tau1).*(1-d1.^(2*initVh))./(1-d1^2);Vp=(d2.^tau2).*(1-d2.^(2*initVp))./(1-d2^2);% standardize the covariance matrices to have determinant = 1 (this avoids% numerical issues of determinants going to infinity or zero)[L,U] = lu(Vh)s =  det(L)        % This is always +1 or -1 diag_u = diag(U)prod_diag_u = prod(diag_u)s_x_prod_diag_u = s*prod_diag_udetVh = det(Vh)detVp = det(Vp);Vh=Vh./detVh^(1/p);Vp=Vp./detVp^(1/q);V=kron(Vp,Vh);invV=V\eye(n);% compute the phylogenetic mean (eq. A.2)U=ones(length(X),1);At=U';Bt=At*invV;fid=fopen('m-out.csv','wt');for i=1:size(Bt,1)  fprintf(fid,'%f ',Bt(i,:));  fprintf(fid,'\n');endfclose(fid);Ct=Bt*UDt=Bt*Xb=(Ct)\(Dt)% compute the phylogenetic MSE (eq. A.2)H=X-b;MSE=(H'*invV*H)/(n-1);