% yCv is the Cv(:) column vector 
% yCv=[Cv11 Cv21 Cv 31 Cv12 Cv22 Cv23 Cv13 Cv23 Cv33]'
% F:= Defromation Gradient
function Cvdot=EvEqn_Cv(t,yCv,x)
global lam

    m1 = x(5); m2 = x(6); a1 = x(7); a2 = x(8); K1 = x(9);
    eta0 = x(10); etai = x(11); beta1 = x(12); beta2 = x(13); K2 = x(14);
    
    
        F = eye(3,3);
        F(3,3) = lam(t);
        F(2,2) = 1/sqrt(F(3,3));  % F11*F22*F33=1
        F(1,1) = 1/sqrt(F(3,3)); 
        det(F);

        
    Cv= reshape(yCv,[3 3]);  % Reshape input yCv into matrix
         
    C=F'*F;

    I1e=trace(C*Cv^-1);
    I1v=trace(Cv);

    I2e=0.5*(I1e^2 - trace((C*Cv^-1)*(C*Cv^-1)));
    J2Neq=((I1e^2)/3-I2e)*(m1*(I1e/3)^(a1-1) + ...
                m2*(I1e/3)^(a2-1))^2;
            

    etao = eta0 ;


    etaK = etai  +( etao - etai + K1*(I1v^beta1 - 3^beta1) )...
                                                /(1 + (K2*J2Neq)^beta2);
                                            
                                            
    eq=(((m1*(I1e/3)^(a1-1) + ...
               m2*(I1e/3)^(a2-1))/etaK)*(C-I1e*Cv/3));
            
    
    Cvdot=eq(:);  % Reshape output as a column vector
end