function [Lep,Nep,Hep] = hyp_stif(y,parms)

% HYP_STIF: computes hypoplastic tensors L and N
%

M=[1,  0,  0,  0,  0,  0
   0,  1,  0,  0,  0,  0
   0,  0,  1,  0,  0,  0
   0,  0,  0,  2,  0,  0
   0,  0,  0,  0,  2,  0
   0,  0,  0,  0,  0,  2];

I=[1,  0,  0,  0,     0,     0
   0,  1,  0,  0,     0,     0
   0,  0,  1,  0,     0,     0 
   0,  0,  0,  0.5,   0,     0
   0,  0,  0,  0,     0.5,   0
   0,  0,  0,  0,     0,   0.5];

kron_delta=[1,  1,  1,  0,  0,  0];

ny = max(size(y));
sig = y(7:12,1);
qint = y(13:ny,1);

e=qint(1);
[p,qdev,z] = inv_s(sig);

T_star=sig-p*kron_delta';
Tstarnorm = sqrt(T_star'*M*T_star);

Mparam     	= parms(4);
lambda       	= parms(5);
kappa     	= parms(6);
Nparam 		= parms(7);
nuparam   	= parms(8);

pcncl=exp((Nparam-log(1+e))/lambda);
psbs=pcncl*(Mparam*Mparam/(Mparam*Mparam+qdev*qdev/(p*p)));
fd=(2*p)/pcncl;
fdsbs=(2*psbs)/pcncl;
fs=3*p/2*(1/kappa+1/lambda)*(1-2*nuparam)/(1+nuparam);

Lep=fs*(I + nuparam/(1-2*nuparam)*kron_delta'*kron_delta);
hypo_Dsom=(3*T_star+kron_delta'*p*(Mparam*Mparam-qdev*qdev/(p*p))/3);
Dsomnorm = sqrt(hypo_Dsom'*M*hypo_Dsom);
hypo_Dsom=hypo_Dsom/Dsomnorm;
	
Amatrix=Lep-sig*kron_delta/lambda;
Nep=Amatrix*hypo_Dsom*fd/fdsbs;

Hep = zeros(6,1);
Hep(1:3) = -(1+e);
