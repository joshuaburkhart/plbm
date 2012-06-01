% IG_AmNat.m

% This is an example program that performs the computations giving the
% strength of phylogenetic signal in host-parasitoid associations.  It was
% used to generate Table 1 in Ives and Godfray "Phylogenetic Analysis of
% Trophic Associations"

% 18 December, 2005

% Required data files
%	hostV.txt - contains the phylogenetic tree of hosts as a covariance
%	matrix

%	paraV.txt - contains the phylogenetic tree of parasitoids as a covariance
%	matrix

%	plantV.txt - contains the phylogenetic tree of plants as a covariance
%	matrix

%	hostparadata.txt - contains the data on host-parasitoid associations as
%	Fik = the number of hosts of species i parasitized by species k in the 
%	first 27 columns and the total number of hosts in column 29 (data
%	from Rott and Godfray, 2000)

% Analyses are performed in 4 ways dictated by the flags:
%	Vhflag=1 => host phylogeny
%	Vhflag=2 => plant phylogeny

%	aflag=1 => use A_ik, the attack rate of parasitoid k on host i (eq 4)
%	aflag=2 => use a_ik, the per capita attack rate (eq 4)

%savefile='bootno01000_3.mat'

%global initVh initVp n p q X tau1 tau2 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize MPI.
MPI_Init;

% Create communicator.
comm = MPI_COMM_WORLD;

% Modify common directory from default for better performance.
% comm = MatMPI_Comm_dir(comm,'/tmp');

% Uncomment if you want to save the messages that were sent.
% comm = MatMPI_Save_messages(comm,1);

% Get size and rank.
comm_size = MPI_Comm_size(comm);
my_rank = MPI_Comm_rank(comm);

% Print rank.
disp(['my_rank: ',num2str(my_rank)]);

% Set who is source and who is destination.
source = 1;
master = 0;

% Create a unique tag id for this message (very important in Matlab MPI!).
tag = 1;


% Seed this rank's rand() function
hostname = getenv('HOSTNAME');
hostnum = sscanf(hostname, '%*[cnu]%i');

% Set up a custom directory for the PCT:
home     = getenv('HOME')
hostpool = sprintf( '%s/.matlab/local_scheduler_data/R2011b/%s' ,home, hostname)
sch = findResource('scheduler','type','local');
set(sch, 'DataLocation', hostpool);

load plantphyno01000_3.txt /ascii
load bactphyno01000_3.txt /ascii
load Ano01000_3.txt /ascii

%matlabpool open 3
matlabpool open 12

%Assign matricies
Vh=plantphyno01000_3(1:end,1:end);
Vp=bactphyno01000_3(1:end,1:end);
A=Ano01000_3(1:end,1:end);


% number of hosts
p=length(Vh);
% number of parasitoids
q=length(Vp);

% Make tips of the phylogenetic trees contemporaneous by extending tips
for i=1:p,
	Vh(i,i)=max(diag(Vh));
end
for i=1:q,
	Vp(i,i)=max(diag(Vp));
end

% scale covariance matrices (this reduces numerical problems caused by
% determinants going to infinity or zero)
Vh=Vh./det(Vh)^(1/p);
Vp=Vp./det(Vp)^(1/q);

X=A;


% Square-root transform data
X=X(:);
X=X.^.5;

n=length(X);

% Compute the Kronecker product of Vp and Vh (eq 1)
V=kron(Vp,Vh);

% Compute the MSE for the case using the "true" phylogenies (d = 1)
U=ones(size(X));
invV=V\eye(n);
aa=(U'*invV*U)\(U'*invV*X);
MSEbase=((X-U*aa)'*invV*(X-U*aa))/(n-1);

% Compute the MSE for the case using the "star" phylogenies (d = 0)
MSEStar=cov(X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EGLS Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let initVh and initVp be the "starter" covariance matrices
initVh=Vh;
initVp=Vp;

% Compute tau
% tau = tau_i + tau_j where tau_i equals the node to tip distance
tau1=diag(initVh)*ones(1,length(initVh)) ...
			+ ones(length(initVh),1)*diag(initVh)'-2*initVh;
tau2=diag(initVp)*ones(1,length(initVp)) ...
			+ ones(length(initVp),1)*diag(initVp)'-2*initVp;


% Numerically find the minimum MSE starting from guess of d in guessd 
options=optimset('MaxFunEvals',10^4);
guessd=[.5 .5];
%[est MSE]=fminsearch('IGfunct',guessd,options);
[est, MSE]=runIGfunctpar(initVh,initVp,n,p,q,X,tau1,tau2, 0, options, guessd);

MSEs=[MSE MSEStar MSEbase]

d1=abs(est(1));
d2=abs(est(2));
d1_d2=[d1 d2]
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		if 1,
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% Create bootstrap data set
			
			% Set the "true" values of X and d
			Xtrue=X;
			dtrue=[d1 d2];

			% Compute covariance matrices with d1 and d2
			Vh=(d1.^tau1).*(1-d1.^(2*initVh))./(1-d1^2);
			Vp=(d2.^tau2).*(1-d2.^(2*initVp))./(1-d2^2);

			% Compute the "true" value of V
			Vh=Vh./det(Vh)^(1/p);
			Vp=Vp./det(Vp)^(1/q);
			Vtrue=kron(Vp,Vh);

			% Compute the "true" value of b
			invV=Vtrue\eye(n);
			btrue=(U'*invV*U)\(U'*invV*Xtrue);

			% Create the matrix T such that T*V*T' = I (identity matrix)
			[L D]=eig(Vtrue);
			iD=diag(diag(D).^(-1/2),0);

			% Construct Y = T*X so that 
			% E{(Y-b)*(Y-b)'} = E{(T*X-b)*(T*X-b)'}
			%				  = T*V*T'
			%				  = I
			
			T=iD*L';
			Y=T*Xtrue;

			% Residuals in orthogonalized space
			resid=Y-btrue;	
			invT=T\eye(n);

			% Number of bootstrap replicates (2000 recommended for 95%
			% confidence intervals)
			reps=60;	
			
      bootlist=zeros(reps,size(dtrue,2)+1);
			parfor i=1:reps,

				% vector of random indices (seed appropriately)
        rng((my_rank+7)*317 + i*67)
				randindex=ceil(n*rand(n,1));	
				
				% create new values of Y selecting residuals at random with
				% replacement
				YY=btrue+resid(randindex);	

				% back-transformed data to create boostrap data set X
				X=invT*YY;	

				% Estimate d for the bootstrap data set
				guessd=dtrue + [0 .1];
				%[est MSE]=fminsearch('IGfunct',guessd,options);
        t = getCurrentTask(); 
        %[est, MSE]=runIGfunctpar(initVh,initVp,n,p,q,X,tau1,tau2, t.ID, options, guessd);
		[est MSE]=mex_gateway(initVh,initVp,n,p,q,X,tau1,tau2,guessd);
				% Collect estimates and MSE into bootlist
				bootlist(i,:) = [abs(est) MSE];		
			end

      if (my_rank > 0 )
        % Send this rank's bootlist.
        MPI_Send( master, tag, comm, bootlist);
      end
      % Destination.
      if (my_rank == master)
        % Receive data.
        for i=1:(comm_size-1)
          dest_boot = MPI_Recv( i, tag, comm );
          bootlist = [bootlist; dest_boot];
          % Check data.
          %if(any((data  - (1:10)) ~= 0))
          %  disp('ERROR: incorrect data sent.');
          %  exit;
          %end
        end
      end

			% Calculate confidence intervals
			%alpha=0.05;
			%conf=[];
			%for i=1:2,
			%	bootlistsorted=sort(bootlist(:,i));
			%	lower_boot=bootlistsorted(floor(0.5*alpha*length(bootlist)));
			%	mean_boot=mean(bootlist(:,i));
			%	upper_boot=bootlistsorted(ceil((1-0.5*alpha)*length(bootlist)));
			%	conf=[conf; lower_boot mean_boot upper_boot];
			%end

			%% Print data
			%lowerconf_mean_upperconf_d=conf

			%min_d=min(bootlist)
			%max_d=max(bootlist)
			%cor_d=corrcoef(bootlist)

			%% Plot bootstrap parameter estimates
			%%figure(2*(aflag-1)+(Vhflag-1)+1)
			%for ii=1:2,
			%	for jj=ii:2,
			%		subplot(2,2,2*(ii-1)+jj)
			%		plot(bootlist(:,ii),bootlist(:,jj),'o',...
			%			dtrue(ii),dtrue(jj),'x')
			%	end
			%end
			%pause(.1)
		end


matlabpool close

%save(savefile);
if (my_rank == master)
  save('pblm_workspace.mat', 'bootlist')
end

% Finalize Matlab MPI.
MPI_Finalize;
disp('SUCCESS');
if (my_rank ~= MatMPI_Host_rank(comm))
  exit;
end
