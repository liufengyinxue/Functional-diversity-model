% last update: 16/7/2015 code comments.
% A Matlab script to calculate the dynamics of a group of species/
% functional groups that compete for water and light.
% This code can be used to produce diversity charts as a function of
% grazing (M1) or precipitation (P), by altering the order of the main
% loops (i.e. hold nm and loop on np for precipitation, hold np and loop on
% nm for grazing).
% Each species is defined by a position on a trait axis, that determines
% both its maximal potential size (K) and its root/shoot ratio (E). 
% Grazing is represented by M1, and the system can be tuned to have
% long/short grazing history with the parameter d0. More explanations and
% the underlying mathematical model are given in the manuscript.
%
% This is supplementary material for the paper: "Linking functional diversity to 
% environmental variability: a mechanistic approach for water limited plant 
% communities".
%
% The code depends on 2 more files that should be in the same folder:
% derivs.m    contains the mathematical formulation of the model equations.
% del2cont.m   return the 2nd derivative of a vector with continuous B.C.
% the results are saved in the file specified in the 'save' command at the
% end of this script, and can be presented using the script
% 'plot_along_p.m'.

% Jonathan Nathan.

clear all; % This is not a function, so we want to make sure the memory is clear.
filename = 'default.mat'; % file name for the data.
% general parameters
time=[0,6000]; % duration of each simulation
nm=2; % run for this # of grazing values
np=100; % run for this # of precipitation values
% tradeoff parameters
emin=0.5; % minimal and maximal values of E, root investment
emax=3.5;
kmin=0.1; % minimal and maximal values of K, shoot investment
kmax=3.5;

% model parameters
% we store all the model parameters in a struct variable called prm.
prm.n=1000; % number of species/functional groups in the simulation
prm.d0=0; % 1=long grazing history (for short d0=0)
prm.kmax=kmax;
alpha=1; % shape of tradeoff curve, as shown in figure 2.
prm.M0=1; % the linear mortality rate
prm.N=2; % the evaporation rate
prm.Lambda_0=0.1; % basic growth rate
prm.h=10; % defines competition for light as in eq. 8 in the paper.
prm.Gamma=1; % water uptake rate by plants
prm.Delta_ek=0.00001; % species mutation rate. 

% create individual traits
xi=linspace(1/prm.n,1-1/prm.n,prm.n); % create a vector for the trait space
prm.E=emin+xi.^alpha*(emax-emin);
prm.K0=kmin+(1-xi).^alpha*((kmax)-kmin);
prm.K=kmin+(1-xi).^alpha*((kmax-prm.d0)-kmin);
prm.d=prm.K0-prm.K; % distance of tradeoff line with grazing history, from grazing line of pure competition.

% initial conditions
prm.x0=[0.01 10]; % initial conditions of biomass and soil-water
x0(1:prm.n)=prm.x0(1);%if we want noise, we add this to I.C.:  +0.1*rand(1,prm.n);
x0(prm.n+1)=prm.x0(2);

% loop parameters
prm.M1_range=linspace(0,4,nm); % range of grazing values to go over
prm.p_range=linspace(400,20,np); % range of precipitation values to go over


if prm.d0 % prevent div.by.0 when there is no grazing history. d0 is changed between 0 and 1. d0 shouldn't be explored along the whole range because of discontinuity around 0.
    prm.beta=kmax/(1.1*prm.d0);% normal value
else 
    prm.beta=0; % value when 0
end

xdat(nm,np,prm.n)=0; % prepare data container for all simulations results.


l=1; % index for lineplot
lineplot(1:nm*np,1:prm.n)=0; % data holder for line plots.

for i=1:nm % run on grazing values
    prm.M1=prm.M1_range(i); % set grazing rate
    for j=1:np % for each grazing value, run on all precipitation values.
        prm.prec=prm.p_range(j); % set prec.
        disp(prm.prec); % show lifesigns
        [t,x]=ode15s(@derivs,time,x0,[],prm); % use a stiff solver, r.h.s file: derivs.m, send parameter struct prm over to derivs function
        x0(1:prm.n)=x(end,1:prm.n); % update I.C. to improve next run's time and result.
        xdat(i,j,:)=(x(end,1:prm.n)); % store solution.
        lineplot(l,:)=(x(end,1:prm.n)); % store lineplot
        btot=sum(x(end,1:prm.n)); % update total biomass value
        rue(i,j)=sum(prm.Lambda_0*(1-(btot-x(end,1:prm.n))/...
            (btot+prm.h)).*x(end,prm.n+1).*x(end,1:prm.n)...
            .*(1+prm.E(1:prm.n).*x(end,1:prm.n)).^2).*(1-x(end,1:prm.n)/prm.K(1:prm.n))/prm.prec; % calculate rain use efficiency and store it.
        l=l+1; % update lineplot index
    end
end
save(filename); % save workspace for further analysis. 
