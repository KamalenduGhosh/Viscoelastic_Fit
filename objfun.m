function [S33,Cv] = objfun(x,data)
global lam
S33 = zeros(length(data),1);
%%

tspan = data(:,1);
Cv0=eye(3,3);      % initialize Cv

opts = odeset('RelTol',1e-8,'AbsTol',1e-12);

[T,Cv] = ode45(@(t,y) EvEqn_Cv(t,y,x),tspan,Cv0(:),opts);  % Pass in column vector initial value
Cv = reshape(Cv.',3,3,[]);             % Reshape the output as a sequence of 3x3 matrices
S = PKStress(T,Cv,x);
S33=squeeze(S(3,3,:));


  
end