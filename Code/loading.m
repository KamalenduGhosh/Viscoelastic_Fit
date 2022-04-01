function lam=loading(t)
global lamdot Fmax Fmax1 Fmax2 Fmax3 Fmax4;
lam_d = lamdot;

tl = (Fmax-1)/lam_d; % time till loading

tl1 = (Fmax1 - 1)/lam_d;
tl2 = (Fmax2 - Fmax1)/lam_d;
tl3 = (Fmax3 - Fmax2)/lam_d;
tl4 = (Fmax4 - Fmax3)/lam_d;
tl5 = (Fmax - Fmax4)/lam_d;
    
% lam=(1+lam_d*t).*(1-heaviside(t-tl))...          % loading branch
%         +(1+2*lam_d*tl-lam_d*t).*heaviside(t-tl);  % unloading branch
% 
lam = (1+lam_d*t).*(1-heaviside(t-tl)) + (1+lam_d*tl)*heaviside(t-tl);   % relaxation

% lam = Fmax1*heaviside(t) + (Fmax2-Fmax1)*heaviside(t-100) + (Fmax4-Fmax3)*heaviside(t-200) +...
%       (Fmax-Fmax4)*heaviside(t-300);

% if t <= 0.1
%     lam = 1 + lam_d*t;
% elseif (t>0.1 & t<= 0.2)
%     lam = Fmax1 + lam_d*(t - 0.1);
% elseif t > (0.2) & (t <= 0.3)
%     lam = Fmax2 + lam_d*(t - 0.2);
% elseif t>(0.3) & t <= (0.4)
%     lam = Fmax3 + lam_d*(t - 0.3);
% elseif t> (0.4) & t <= (0.5)
%     lam = Fmax4 + lam_d*(t - 0.4);
% else 
%     lam = Fmax;
    


end