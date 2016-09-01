%% Parameters for the particle filter
Np = 2e3; %number of particles in filter

%OBSERVATIONAL ERROR COVARIANCES
u_sigma=2e-1; %observational error variance
v_sigma = u_sigma;
t_sigma = .5;
%OBSERVATION COVARIANCE MATRIX SET AT END OF FILE

%noise in parameter values after resampling
lw = 0.01;
aw = 0.01;
tw = 0.1;

resamp_thresh=0.4; %threshold for resampling
 

%% Generate particles for PF
% PARTICLE IS [L_INF, ALPHA, THETA_INF]
% L_INF ~ 10-20
% ALPHA ~ 4-5
% THETA_INF ~ 290

min_guess=[8;3;270]; max_guess=[22;6;320]; %in column

dim_param = length(min_guess);
particle = samp_param(min_guess,max_guess,Np);
W = ones(Np,1)/Np;  %uniform weight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate observations, initial conditions etc for PDE model

%NOTE: MOST PHYSICAL CONSTANTS DEFINED IN INTEGRATION CODE

nn=15; %which day to start

%timestep
dt=10; %measured in seconds
t = 0;          %min value of t
hr = 3600; %seconds in one hour
min = 60; %seconds in one minute
t_int = 10*min;     % observation interval
nPF = 54; %number of steps in the particle filter; the total time integrated over the PF will be nPF * t_int
%MAXIMUM VALUE OF nPF IS 54, UNTIL WE STITCH TOGETHER DIFFERENT NIGHTS


x = sort([ 20:10:500 10 2 .1 ])'; %really fine grid
%x = sort([ 500 200 140 100 40 20 10 2 .1 ])'; %minimum grid
%x = sort([ 200:60:500 140 100 40 20 10 2 .1 ])'; %finer grid at large x

%x = refine(x,1);
L=length(x);
N = L-1;
%ghost nodes
xl=(2*x(1))-x(2);
xr=x(L)+(x(L)-x(L-1));
x=[xl; x; xr];

xobs=[ 500 200 140 100 40 20 10 2 .1 ]; %observation locations, plus boundaries

% INDICES FOR DATA POINTS
data_pts =[2 10 20 40 100 140 200]; %2 not included for u,v
It = zeros(length(data_pts),1); %indices for theta data
for j=1:length(data_pts)
    It(j) = find(x==data_pts(j)); %find the indices in x that correspond to observation locations for theta
end
Iuv = It(2:end); %observation locations for u, v




nsteps = t_int/dt;

%define variables
xMin = 0.1;   xMax=500;       %max and min point on the main grid

% %generate main grid (nonuniform)
% xc1 = logspace(log10(xMin),log10(2), 4);
% xc2 = logspace(log10(2),log10(10),4);
% xc3 = logspace(log10(10),log10(80),4);
% xc4 = 140:60:xMax;
% x=[xc1 xc2(2:end) xc3(2:end) xc4]'; %The 2:end entries delete repeated entries
% %refine grid here
% x=refine(x,2);





% STAGGERED GRID DEFINED IN INTEGRATION CODE


%note that 
%1. 'x' and 'xs' are nonuniform grids
%2. 'x' and 'xs' both include ghost nodes




% geowind data

datatemp=importdata('data_for_kansas_geowind.mat');

%%%%%%%%%%%%%%
start_3=4365; %start location for third season
sta=start_3+24*(nn-1); %geowind data is hourly

geowind=datatemp(4,:);
direc=datatemp(5,:)*pi/180;

ug=-geowind.*sin(direc);
UgM=ug(sta:sta+9); %includes 10 data points due to interpolation in the loop. The tenth data point does not affect the integration.
vg=-geowind.*cos(direc);
VgM=vg(sta:sta+9); %includes 10 data points due to interpolation in the loop


% initial conditions
load('data_for_kansas_fluxes.mat','SH_ssn') %load only SH_ssn because otherwise variables like dt are overwritten
load('data_for_kansas_tower.mat')


%season 3
u3 = u_long_ssn{3};
v3 = v_long_ssn{3};
t3 = theta_long_ssn{3};
un3 = u3(night_pts_ssn{3},:);
vn3 = v3(night_pts_ssn{3},:);
tn3 = t3(night_pts_ssn{3},:);
%refine vectors to only the current night
in_p = 1 + (nn-1)*54; %location of initial conditions in matrices; nights are 54 data points, days are 90
un3 = un3(in_p:in_p+54,1:end-1); % 
vn3 = vn3(in_p:in_p+54,1:end-1);
tn3 = tn3(in_p:in_p+54,:);

Hs = SH_ssn{3}(night_pts_ssn{3})/1200; %scaled to be RHS of boundary condition
Hs = Hs(in_p:in_p+54); %truncate Hs and return values for one night only

% solve Km * theta_z = Hs(1)

%
uobs = [UgM(1),un3(1,1:end),NaN,0]; %initially observed values: geowind at boundary, data in between, no observation at z=2m, no wind at z=0.1m
vobs = [VgM(1),vn3(1,1:end),NaN,0];


%set initial conditions
u = zeros(1,length(x));
v=u;
theta = zeros(Np,length(x));
%M(1,:) = interp1([xobs xl],[uobs 0],x,'spline'); %sets u to 0 at ghost node 
u(1,:) = interp1(xobs,uobs,x,'spline');

u = repmat(u,Np,1); %Now u(n,:) is the forecast for u for the n-th particle

%M(2,:) = interp1([xobs xl],[vobs 0],x,'spline'); %sets v to 0 at ghost node
%forecast = repmat(M,N,1);
v(1,:) = interp1(xobs,vobs,x,'spline');

v = repmat(v,Np,1);

% INITIAL CONDITIONS FOR THETA DEPEND ON THE VALUE OF THETA_INF
tntemp = tn3(1,:);
parfor n=1:Np
    theta_inf = particle(3,n);
    tobs = [theta_inf,tntemp];
    theta(n,:) = interp1([ 500 200 140 100 40 20 10 2],tobs,x,'spline');
end
clear theta_inf tntemp %do not retain misleading variables





len = length(data_pts)-1; %number of u data points
R=eye(3*len+1); %observational error covariance matrix. Note independence assumption in using a diagonal matrix.
for j = 1:len
    R(j,j) = u_sigma;
    R(len+j,len+j) = v_sigma;
    R(2*len+j,2*len+j) = t_sigma;
end
R(end,end) = t_sigma; %there is always one more data point for theta than u,v





   % UgM(1)
   % VgM(1)
   
   
%N = # of main nodes -1


%  u(2)
%  u(N+2)
%  theta(2)
%  theta(N+2)
%aviobj = avifile('scm.avi','compression','None');



%TEST FOR ONE PARTICLE
% n = 1+round((N-1)*rand); %choose a random particle
% nfor = 1+3*(n-1):3*n; %range for current forecast
% param = particle(:,n);
% tic
% [u,v,theta] = int_pde_n([u;v;theta(n,:)],x,nsteps,dt,UgM,VgM,N,Hs,param);
% toc
% figure
% hold on
% plot(x(2:N+2),u(2:N+2),'rx-');
% plot(x(2:N+2),v(2:N+2),'rx-');
% plot(x(2:N+2),theta(2:N+2),'rx-')