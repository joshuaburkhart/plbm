% Computes the MSE for the model (eq 2)%	A = Xb + ?function MSE=IGfunct(parameter)global initVh initVp n p q X tau1 tau2d1=parameter(1);d2=parameter(2);invV=make_Inverse(d1 d2);MSE=make_MSE(invV);