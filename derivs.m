function dydx=derivs(t,x,prm) % rhs for the model equations.
dydx(1:prm.n+1)=0; % initialize result variable
uptake(1:prm.n+1)=0; % temporary holder for integration over uptake of all species
% calculate general terms in the loop
d2x(1:prm.n)=del2cont(x(1:prm.n),prm.n); % mutation: discrete second derivative in trait space.
btot=sum(x(1:prm.n)); % usefull value for later calculations.
fw=x(end); % water level at time of calculation.

for i=1:prm.n % run over all species
    m=prm.M0*x(i)+prm.M1*(1-prm.beta*prm.d(i)/prm.kmax)*x(i).^2; % calculate mortality term, include grazing history and small correction "beta" 
    Lambda=prm.Lambda_0*(1-(btot-x(i))/(btot+prm.h));
%     Lambda=prm.Lambda_0;
    uptake(i)=fw*(1+prm.E(i)*x(i)).^2;
    comp=x(i)/prm.K(i);
    dydx(i) = Lambda*x(i)*uptake(i).*(1-comp)-m + prm.Delta_ek*d2x(i); % species equation, 1a in the paper
end
dydx(prm.n+1) = prm.prec - prm.N*x(end) - sum(prm.Gamma*x(1:prm.n).*uptake(1:prm.n)') ; % water equation, 1b in the paper
dydx=dydx'; % return a vector of results.
