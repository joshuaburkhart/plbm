% Computes the MSE for the model (eq 2)%	A = Xb + ?function MSE=IGfunct(parameter)global initVh initVp n p q X tau1 tau2invV=make_Inverse(parameter);MSE=make_MSE(invV);