% yCv is the Cv(:) column vector 
% yCv=[Cv11 Cv21 Cv 31 Cv12 Cv22 Cv23 Cv13 Cv23 Cv33]'
% F:= Defromation Gradient
function [J2 I1v etaK] = Findeta(t,yCv,x,type)
    global nu1 nu2 b1 b2 eta0 conc lam
    
    if type==0  % Uniaxial
        F = eye(3,3);
        F(3,3) = lam(t);
        F(2,2) = 1/F(3,3);  % F11*F22*F33=1
        F(1,1) = 1; 
        det(F);
    elseif type ==1 %Simple Shear
        F = eye(3,3);
        F(1,2) = lam(t)-1;
        det(F);
    end
        
    Cv= reshape(yCv,[3 3]);  % Reshape input yCv into matrix
         
    C=F'*F;

    I1e=trace(C*Cv^-1);
    I1v=trace(Cv);

    I2e=0.5*(I1e^2 - trace((C*Cv^-1)*(C*Cv^-1)));
    J2Neq=((I1e^2)/3-I2e)*(nu1*(I1e/3)^(b1-1) + ...
                nu2*(I1e/3)^(b2-1))^2;
            
    J2 = 0.5*(nu1^2)*(trace(Cv^-1*C*Cv^-1*C) - ...
                      trace((C*Cv^-1))^2/3);

    etaR = eta0*(1-conc);
    etai = etaR/(1 + (x(1)*(I1v)/nu1));
    etaK = etai + (eta0 - etai)/(1 + (x(1)*J2/nu1^2)^x(1));

    eq=((nu1*(I1e/3)^(b1-1) + ...
                nu2*(I1e/3)^(b2-1))/etaK)*(C-I1e*Cv/3);
    
    Cvdot=eq(:);  % Reshape output as a column vector
end