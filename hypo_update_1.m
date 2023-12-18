function y_kp1 = hypo_update_1(y,V,parms,k)

% HYPO_UPDATE: integrate evolution equations of hypoplasticity model
%

tiny = 1.0e-10;

M = [1,  0,  0,    0,    0,    0
     0,  1,  0,    0,    0,    0
     0,  0,  1,    0,    0,    0
     0,  0,  0,    2,    0,    0
     0,  0,  0,    0,    2,    0
     0,  0,  0,    0,    0,    2];

%% recover input data

max_ksub = parms(1);
err_tol  = parms(2);

ny = max(size(y));

%% initialization

T_j=0.0;
DT_j=1.0;
y_j=y;

%% start adaptive substepping for the integration of the elastoplastic eqs.

ksub=0;

while T_j<1.0		

	ksub=ksub+1;

% stop if too many substeps

    if ksub > max_ksub
        disp ('--- ERROR: substep number exceeding maximum in trial state evaluation ---')
		error('    execution stopped in HYPO_UPDATE_1')
    end
    
% build Runge-Kutta functions

	[E,S] =  constraints(k,y_j);
    KRK_1=f_hyp(y_j,S,E,V,parms);
    
	y_2=y_j+(DT_j/2.0)*KRK_1;

	[E,S] = constraints(k,y_2);
	KRK_2=f_hyp(y_2,S,E,V,parms);
    
	y_3=y_j-DT_j*KRK_1+2.0*DT_j*KRK_2;

    [E,S] = constraints(k,y_3);
	KRK_3=f_hyp(y_3,S,E,V,parms);

% form approximate solution of 2nd and 3rd order
    
	y_hat=y_j+DT_j*KRK_2;
	y_til=y_j+DT_j*((1.0/6.0)*KRK_1+(2.0/3.0)*KRK_2+(1.0/6.0)*KRK_3);

% local error estimate

    sig_til  = y_til(7:12,1);
    qint_til = y_til(13:ny,1);
    
    delta_sig  = y_til(7:12,1)-y_hat(7:12,1);
    delta_qint = y_til(13:ny,1)-y_hat(13:ny,1);
    
    norm_sig  = sqrt(sig_til'*M*sig_til);
    norm_qint = sqrt(qint_til'*qint_til);
      
    if norm_sig < tiny      
        norm_sig=tiny;
    end
    
    if norm_qint < tiny
        norm_qint=tiny;
    end

    Res = zeros(ny,1);
    Res(7:12,1)=(1.0/norm_sig)*delta_sig;
    Res(13:ny,1)=(1.0/norm_qint)*delta_qint;

% residual vector and normalized step size
    
    norm_Res=sqrt(Res'*Res);
    
    if norm_Res < tiny;
        norm_Res = tiny;
    end

	NSS = 0.9*DT_j*(err_tol/norm_Res)^(1.0/3.0);

	if norm_Res<err_tol 				

% substep is accepted, update y and T and estimate new substep size DT

		label=' -- hypoplastic integration -- accepted; T = ';
				
		y_j=y_til;
					
		T_j=T_j+DT_j;
		DT_j=min(4.0*DT_j,NSS);
		DT_j=min((1.0-T_j),DT_j);
%          disp([label,num2str(T_j,4),' - DT = ',num2str(DT_j,4)])

	else

% substep is rejected, recompute with new (smaller) substep size DT

		label=' -- hypoplastic integration -- rejected; T = ';
								
		DT_j=max(0.25*DT_j,NSS);
%          disp([label,num2str(T_j,'single'),' - DT = ',num2str(DT_j,'single')])
					
	end						

% bottom of while loop

end 

% update state variables at the end of substepping

y_kp1 = y_j;

