% INPUT_DATA: input material constants and store them in vector parms
%             set loading path info
%
% nspb          = no. of stress path branches;
% icode(1:nspb) = vector of path type codes (one for each branch)
% nstep(1:nspb) = vector of step numbers (one for each branch)
% DX(1:nspb)    = vector of total "load" increment for each branch
%
% loading path codes
%
% 1. strain controlled undrained TX compression (axis direction: x_3)
% 2. strain controlled ED compression (axis direction: x_3)
% 3. stress controlled drained TX compression (axis direction: x_3)
% 4. mixed control drained TX compression
% 5. mixed control ED compression
% 6. mixed control plane strain compression
% 7. strain controlled undrained simple shear
% 8. stress (mixed) controlled simple shear

%% material parameters

Mparam     	= 1;
lambda       	= 0.1;
kappa     	= 0.01;
Nparam 		= 1;
nuparam   	= 0.2;

max_ksub 	= 1000;       % max no. of substeps allowed
err_tol  	= 1.0e-4;     % normalized error tolerance for adaptive integration  
tol_mnriter	= 1.0e-10;     % tolerance on Modified Newton Raphson estimation of deps

parms=[max_ksub;
       err_tol;
       tol_mnriter;
       Mparam;
       lambda;
       kappa;
       Nparam;
       nuparam];

%% type and characteristics of loading path

nspb=1;

icode=[4];
nstep=[100];
DX=[0.15];

path_info = [icode;
             nstep
             DX];
