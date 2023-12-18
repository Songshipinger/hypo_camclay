function F = f_hyp(y,S,E,V,parms)

% F_EP: computes RHS of the evolution ODE for plastic processes
%

M=[1,  0,  0,  0,  0,  0
   0,  1,  0,  0,  0,  0
   0,  0,  1,  0,  0,  0
   0,  0,  0,  2,  0,  0
   0,  0,  0,  0,  2,  0
   0,  0,  0,  0,  0,  2];

M2=[1,  0,  0,  0,    0,    0
   0,  1,  0,  0,    0,    0
   0,  0,  1,  0,    0,    0 
   0,  0,  0,  0.5,  0,    0
   0,  0,  0,  0,    0.5,  0
   0,  0,  0,  0,    0,    0.5];


ny = max(size(y));
F = zeros(ny,1);

[Lep,Nep,Hep] = hyp_stif(y,parms);

%% compute A matrix and solve for deps

Laux=2*Lep;
A    = S*Laux+E;
depsMNR = A\V;
depsinit = depsMNR;
sigpresc = Laux*depsMNR;

%% Modified Newton-Raphson interations for correction of deps
stresscontrol=zeros(6,1);
for i = 1:6
  for j = 1:6
    if  (S(i,j) > 0)
      stresscontrol(j)=1;
    end
  end
end

isstresscont=sqrt(stresscontrol'*M*stresscontrol);

if isstresscont > 0
  tol_mnriter  = parms(3);
  itererror=1000;
  
  while itererror > tol_mnriter
    
        sigest= Lep*depsMNR - Nep*sqrt(depsMNR'*M2*depsMNR);
        unbalsig=sigpresc-sigest;
 
        for i = 1:6
	  if  (stresscontrol(i) == 0)
	    depsMNR(i)=depsinit(i);
	    unbalsig(i)=0;
	  end
        end

        depsMNR=depsMNR + Laux\unbalsig;
        itererror=sqrt(unbalsig'*M*unbalsig);
        
  end
end

%% F vector

F(1:6,1)   = depsMNR;
F(7:12,1)  = Laux*depsMNR - Nep*sqrt(depsMNR'*M2*depsMNR);
de=Hep'*depsMNR;
F(13:ny,1) = de;
