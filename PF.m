%This file runs particle filter on a supplied model and parameters
clear
clc
PF_param; %runs m-file to generate model and model parameters



% Loop over time
for tau=1:nPF
   tau %tau is an index; (tau-1) * nsteps * dt is the initial time for each integration
    obs = [un3(tau,:) vn3(tau,:) tn3(tau,:)]; %get observations for current time
    % Get observations
    
    %particle filter
    parfor n=1:Np
        nfor = 1+3*(n-1):3*n; %range for current forecast
        param = particle(:,n);
        [ut,vt,thetat] = int_pde_n([u(n,:);v(n,:);theta(n,:)],x,nsteps,dt,UgM,VgM,N,Hs,param);

        u(n,:) = ut;
        v(n,:) = vt;
        theta(n,:) = thetat;
        xend = [ut(Iuv)' vt(Iuv)' thetat(It)'];
        W(n) = W(n)*exp(-1/2*(obs-xend)/R*(obs-xend)'); %note Gaussian assumption
    end
    W(isnan(W))=0;
    W=W/sum(W);
    
    
    %test for filter collapse
    if isnan(W(1)) 
        error('Particle filter collapsed at time step t = %f',tau)
    %test for resampling
    else if 1/sum(W.^2)/Np < resamp_thresh 
            sampIndex = ResampSimp(W,Np);
            particle = particle(:,sampIndex)+[lw*randn(1,Np);aw*randn(1,Np);tw*randn(1,Np)];
            u = u(sampIndex,:);
            v = v(sampIndex,:);
            theta = theta(sampIndex,:);
            W=1/Np*ones(Np,1);
        end
    end
    
end

%% Graphs and Plotting

 figure
 hold on
 histogram(particle(1,:),50)%,'Normalization','Probability'),'Normalization','Probability')
 xlabel('l_\infty','FontSize',15)
 ylabel('Particle count','FontSize',15)
 
 
 figure
 histogram(particle(2,:),50)%,'Normalization','Probability'),'Normalization','Probability')
 xlabel('\alpha','FontSize',15)
  ylabel('Particle count','FontSize',15)

  figure
 histogram(particle(3,:),50)%,'Normalization','Probability')
 xlabel('\theta_\infty','FontSize',15)
  ylabel('Particle count','FontSize',15)


figure
plot3(particle(1,:),particle(2,:),particle(3,:),'o','MarkerSize',5)
xlabel('l_\infty','FontSize',15)
ylabel('\alpha','FontSize',15)
zlabel('\theta_\infty','FontSize',15)

