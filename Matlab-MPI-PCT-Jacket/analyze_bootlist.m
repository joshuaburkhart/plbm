load pblm_workspace
dtrue = [0.2106    0.8695];


			% Calculate confidence intervals
			alpha=0.05;
			conf=[];
			for i=1:2,
				bootlistsorted=sort(bootlist(:,i));
				lower_boot=bootlistsorted(floor(0.5*alpha*length(bootlist)));
				mean_boot=mean(bootlist(:,i));
				upper_boot=bootlistsorted(ceil((1-0.5*alpha)*length(bootlist)));
				conf=[conf; lower_boot mean_boot upper_boot];
			end

			% Print data
			lowerconf_mean_upperconf_d=conf

			min_d=min(bootlist)
			max_d=max(bootlist)
			cor_d=corrcoef(bootlist)

			% Plot bootstrap parameter estimates
			%figure(2*(aflag-1)+(Vhflag-1)+1)
			for ii=1:2,
				for jj=ii:2,
					subplot(2,2,2*(ii-1)+jj)
					plot(bootlist(:,ii),bootlist(:,jj),'o',...
						dtrue(ii),dtrue(jj),'x')
				end
			end
			pause(.1)

