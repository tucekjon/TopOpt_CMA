function [xmma, xold1, xold2, low, upp, AsymInc, AsymDecr, change] = ...
          mma_unconstrained(dfdx, xval, xold1, xold2, ...
                      low, upp, iter, AsymInc, AsymDecr, xmin, xmax)  


move = 0.2;
Osc = 0.2; 
AsymInit = 0.5;

xmin = max(xmin,xval-move); 
xmax = min(xmax,xval+move);

% Compute asymptotes L and U:
AsymInc = min(1 + Osc, AsymInc); 
AsymDecr = max(1 - 2*Osc, AsymDecr);

if iter<=2
  low = xval - AsymInit*(xmax-xmin);  
  upp = xval + AsymInit*(xmax-xmin);    
else
  sgn = (xval-xold1).*(xold1-xold2);
  s = ones(size(xval)); 
  s(sgn>0) = AsymInc; 
  s(sgn<0) = AsymDecr;

  low = xval - s.*(xold1 - low); 
  upp = xval + s.*(upp - xold1);
end

% Compute bounds alpha and beta
alpha = 0.9*low + 0.1*xval;   
beta = 0.9*upp + 0.1*xval;

alpha = max(xmin, alpha); 
beta = min(xmax, beta);

% Solve unconstrained subproblem
feps = 0.000001; 
p = (upp-xval).^2.*(max(dfdx,0)+0.001*abs(dfdx)+feps./(upp-low)); 
q = (xval-low).^2.*(-min(dfdx,0)+0.001*abs(dfdx)+feps./(upp-low));
zCnd = (low.*p - upp.*q + (upp - low).*sqrt(p.*q))./(p - q);
xmma = max(alpha,min(beta,zCnd));

xold2 = xold1; 
xold1 = xval;
% change = sum(abs(xmma-xval))/length(xval);
change = max(abs(xmma-xval));