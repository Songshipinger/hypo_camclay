% INIT_STATE: initialize state variables
%
% eps0  = initial strain state
% sig0  = initial stress state
% qint0 = initial value of internal variables (stored in a 1d array)

%% preallocate arrays

eps0 = zeros(6,1);
sig0 = zeros(6,1);

%% input initial strain

% eps0 = [0.0; 0.0; 0.001; 0.0; 0.0; 0.0]   % uncomment this line if eps0 =/= 0

%% input initial stress

sig0 = [100.0; 100.0; 100.0; 0.0; 0.0; 0.0];
e0 = 0.46;

%% input initial internal variables

qint0 = [e0];

%% collect all state variables in vector y0

y0 = [eps0',sig0',qint0']';
