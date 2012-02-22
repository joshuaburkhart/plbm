% IG_AmNat.m% This is an example program that performs the computations giving the% strength of phylogenetic signal in host-parasitoid associations.  It was% used to generate Table 1 in Ives and Godfray "Phylogenetic Analysis of% Trophic Associations"% 18 December, 2005% Required data files%	hostV.txt - contains the phylogenetic tree of hosts as a covariance%	matrix%	paraV.txt - contains the phylogenetic tree of parasitoids as a covariance%	matrix%	plantV.txt - contains the phylogenetic tree of plants as a covariance%	matrix%	hostparadata.txt - contains the data on host-parasitoid associations as%	Fik = the number of hosts of species i parasitized by species k in the %	first 27 columns and the total number of hosts in column 29 (data%	from Rott and Godfray, 2000)% Analyses are performed in 4 ways dictated by the flags:%	Vhflag=1 => host phylogeny%	Vhflag=2 => plant phylogeny%	aflag=1 => use A_ik, the attack rate of parasitoid k on host i (eq 4)%	aflag=2 => use a_ik, the per capita attack rate (eq 4)global initVh initVp n p q X tau1 tau2load hostV.txt /asciiload paraV.txt /asciiload plantV.txt /asciiload hostparadata.txt /asciifor aflag=1:2,	for Vhflag=1:2,		flags=[aflag,Vhflag]		% Assign host phylogeny		if Vhflag==1,			Vh=hostV(2:end,2:end);		end		if Vhflag==2,			Vh=plantV(2:end,2:end);		end		% Assign parasitoid phylogeny		Vp=paraV(2:end,2:end);		% number of hosts		p=length(Vh);		% number of parasitoids		q=length(Vp);		% Make tips of the phylogenetic trees contemporaneous by extending tips		for i=1:p,			Vh(i,i)=max(diag(Vh));		end		for i=1:q,			Vp(i,i)=max(diag(Vp));		end		% scale covariance matrices (this reduces numerical problems caused by		% determinants going to infinity or zero)		Vh=Vh./det(Vh)^(1/p);		Vp=Vp./det(Vp)^(1/q);		% Assign data		F=hostparadata(:,1:27);		hostden=hostparadata(:,29);		% Equation (4)		H=hostden*ones(1,27);		A=-log(1-F./H);		% Use A_ik		if aflag==1,			X=A;		end		% Use a_ik by standardizing across parasitoid densities		if aflag==2,			meanA=[];			for i=1:q,				meanA=[meanA mean(A(find(A(:,i)>0),i))];				%meanA=[meanA mean(A(:,i))];			end			meanst=ones(p,1)*meanA;			X=A./meanst;		end		% Square-root transform data		X=X(:);		X=X.^.5;		n=length(X);		% Compute the Kronecker product of Vp and Vh (eq 1)		V=kron(Vp,Vh);		% Compute the MSE for the case using the "true" phylogenies (d = 1)		U=ones(size(X));				%invV=V\eye(n); %this takes a long time		Vg=gpuArray(V); %push V to gpu memory		invVg=Vg\parallel.gpu.GPUArray.eye(n); %try this instead		invV=gather(invVg); %push back to cpu				aa=(U'*invV*U)\(U'*invV*X);		MSEbase=((X-U*aa)'*invV*(X-U*aa))/(n-1);		% Compute the MSE for the case using the "star" phylogenies (d = 0)		MSEStar=cov(X);		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		% EGLS Estimation		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		% Let initVh and initVp be the "starter" covariance matrices		initVh=Vh;		initVp=Vp;		% Compute tau		% tau = tau_i + tau_j where tau_i equals the node to tip distance		tau1=diag(initVh)*ones(1,length(initVh)) ...			+ ones(length(initVh),1)*diag(initVh)'-2*initVh;		tau2=diag(initVp)*ones(1,length(initVp)) ...			+ ones(length(initVp),1)*diag(initVp)'-2*initVp;		% Numerically find the minimum MSE starting from guess of d in guessd 		options=optimset('MaxFunEvals',10^4);		guessd=[.5 .5];		[est MSE]=fminsearch('IGfunct',guessd,options);		MSEs=[MSE MSEStar MSEbase]		d1=abs(est(1));		d2=abs(est(2));		d1_d2=[d1 d2]			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		% Bootstrap		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		if 1,						%%%%%%%%%%%%%%%%%%%%%%%%%%%			% Create bootstrap data set						% Set the "true" values of X and d			Xtrue=X;			dtrue=[d1 d2];			% Compute covariance matrices with d1 and d2			Vh=(d1.^tau1).*(1-d1.^(2*initVh))./(1-d1^2);			Vp=(d2.^tau2).*(1-d2.^(2*initVp))./(1-d2^2);			% Compute the "true" value of V			Vh=Vh./det(Vh)^(1/p);			Vp=Vp./det(Vp)^(1/q);			Vtrue=kron(Vp,Vh);			% Compute the "true" value of b			invV=Vtrue\eye(n);			btrue=(U'*invV*U)\(U'*invV*Xtrue);			% Create the matrix T such that T*V*T' = I (identity matrix)			[L D]=eig(Vtrue);			iD=diag(diag(D).^(-1/2),0);			% Construct Y = T*X so that 			% E{(Y-b)*(Y-b)'} = E{(T*X-b)*(T*X-b)'}			%				  = T*V*T'			%				  = I						T=iD*L';			Y=T*Xtrue;			% Residuals in orthogonalized space			resid=Y-btrue;				invT=T\eye(n);			% Number of bootstrap replicates (2000 recommended for 95%			% confidence intervals)			reps=2000;							%bootlist=[];			bootlist=zeros(reps,3);			%for i=1:reps,			parfor i = 1:reps    % what the hell.. give it a shot				% vector of random indices				randindex=ceil(n*rand(n,1));									% create new values of Y selecting residuals at random with				% replacement				YY=btrue+resid(randindex);					% back-transformed data to create boostrap data set X				X=invT*YY;					% Estimate d for the bootstrap data set				guessd=dtrue + [0 .1]				guessd				%%%%%%%%%				exit;				%%%%%%%%%				[est MSE]=fminsearch('IGfunct',guessd,options);			%%%%%%%%%%%%%%%%%%%%%%%%%			%abs(est)			%MSE			%class(est)			%class(abs(est))			%class(MSE)			%size(abs(est))			%size(MSE)			%%%%%%%%%%%%%%%%%%%%%%%%%				% Collect estimates and MSE into bootlist				bootlist(i,:)=[abs(est) MSE];								%%%%%%%%%%%%%%%%%%%%				%bootlist				%%%%%%%%%%%%%%%%%%%%							end			% Calculate confidence intervals			alpha=0.05;			conf=[];			for i=1:2,				bootlistsorted=sort(bootlist(:,i));				lower_boot=bootlistsorted(floor(0.5*alpha*length(bootlist)));				mean_boot=mean(bootlist(:,i));				upper_boot=bootlistsorted(ceil((1-0.5*alpha)*length(bootlist)));				conf=[conf; lower_boot mean_boot upper_boot];			end			% Print data			lowerconf_mean_upperconf_d=conf			min_d=min(bootlist)			max_d=max(bootlist)			cor_d=corrcoef(bootlist)			% Plot bootstrap parameter estimates			figure(2*(aflag-1)+(Vhflag-1)+1)			for ii=1:2,				for jj=ii:2,					subplot(2,2,2*(ii-1)+jj)					plot(bootlist(:,ii),bootlist(:,jj),'o',...						dtrue(ii),dtrue(jj),'x')				end			end			pause(.1)		end	endend