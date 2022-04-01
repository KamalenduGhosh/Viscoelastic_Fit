clc; clear;
global lam lamdot Fmax 
lamdot = 0.05;
Fmax = 3;

formatSpec = ['%f,%f'];
fname = strcat('LoadingUnloading_F',num2str(Fmax),'_l',num2str(lamdot),'.txt');
fileID = fopen(fname,'r');
Data = fscanf(fileID,formatSpec,[2,Inf])';
fclose(fileID);

figure(1)
plot(Data(:,1), Data(:,2),'.')
title({'$\textrm{Data for VHB }4910,|\dot{\lambda}|=0.05 s^-1 $'},'Interpreter','latex','FontSize',15)
ylabel({'$S_{33}(\textrm{kPa})$'},'Interpreter','latex','FontSize',15)
xlabel({'$F_{33}=\lambda$'},'Interpreter','latex','FontSize',15)
%% Convert stretch data into time history data
Data = [zeros(length(Data(:,1)),1) Data(:,1) Data(:,2)];

ind1 = find(Data(:,2) == max(Data(:,2)));
ind2 = find(Data(:,2) == min(Data(:,2)));
for i = 1:length(Data(:,2))
    if i <= ind1
        Data(i,1) = (Data(i,2) - 1)/lamdot;
    else
        tload = (Fmax - 1)/lamdot;
        Data(i,1) = (1  + 2*lamdot*tload - Data(i,2))/lamdot;
    end
end
% % 
lam = @(t) (1+lamdot*t).*(1-heaviside(t-tload))+...
        ((1+2*lamdot*tload-lamdot*t).*heaviside(t-tload));

figure (2)
plot(Data(:,1) , lam(Data(:,1)))
title({'$\textrm{Loading History } $'},'Interpreter','latex','FontSize',15)
ylabel({'$\lambda$'},'Interpreter','latex','FontSize',15)
xlabel({'$t(\textrm{s})$'},'Interpreter','latex','FontSize',15)
%% Solve Evolution Equation using RK5

lb = [];
ub = [];
A = [];
b = [];
Aeq = [];
beq = [];
opt = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Display','iter',...
    'FunctionTolerance', 1e-12,'StepTolerance', 1e-12,'MaxFunctionEvaluations',1,'MaxIterations',0);


% initial guess 
x0 = [13.54 1.08 1.00 -2.474 5.42 20.78 -10 1.948 3507 7014 0.1 1.852 0.26 1];
%%

[x,resnorm,residual,exitflag,output] = lsqnonlin(@(x) norm(objfun(x,Data(1:end,:)) - Data(1:end,3)),x0, lb, ub,opt);

Sfitted = objfun(x,Data(1:end,:));
figure(7)
plot(Data(1:end,2),Data(1:end,3),'b.','MarkerSize',12)
hold on;
plot(Data(1:end,2),Sfitted,'r-','LineWidth',2)
grid on;

title({'$\textrm{VHB }4910,|\dot{\lambda}|=0.05 s^-1 $'},'Interpreter','latex','FontSize',15)
ylabel({'$S_{33}(\textrm{kPa})$'},'Interpreter','latex','FontSize',15)
xlabel({'$F_{33}=\lambda$'},'Interpreter','latex','FontSize',15)
legend({'$\textrm{Data}$',...
     '$\textrm{Model } $'},...
     'Interpreter','latex','FontSize',10)
hold off
