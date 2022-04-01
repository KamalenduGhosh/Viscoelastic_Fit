function S=PKStress(t,Cv,x)
global lam
    mu1 = x(1); mu2 = x(2); alpha1 = x(3); alpha2 = x(4);
    m1 = x(5); m2 = x(6); a1 = x(7); a2 = x(8);
    
    %%
    S=zeros(3,3,length(t));
    F=zeros(3,3,length(t));
    
 
    F(3,3,:) = lam(t);
    F(2,2,:) = 1/sqrt(F(3,3,:));
    F(1,1,:) = 1/sqrt(F(3,3,:));
    
        
    
    for i=1:length(t)
        C=F(:,:,i)'*F(:,:,i);
        I1=trace(C);
        I1e=trace(C*Cv(:,:,i)^-1);
        
        S1=(mu1*(I1/3)^(alpha1-1) + mu2*(I1/3)^(alpha2-1))*F(:,:,i)+...
            (m1*(I1e/3)^(a1-1) + m2*(I1e/3)^(a2-1))*F(:,:,i)*Cv(:,:,i)^(-1);
        
        p=S1(2,2)*F(2,2,i);  % Finding p
        
        S(:,:,i)=S1-p*(F(:,:,i)^-1)'; 
           
    end   
end