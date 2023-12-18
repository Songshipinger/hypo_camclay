% DRIVER_EP - Octave version
%
% Constitutive driver for integrating elastoplastic constitutive
% equations using the algorithm proposed by
% Bardet & Choucair (1991), IJNAMG
%
% written by C. Tamagnini - Jan. 2012
%

clear;
close all;
clc;

disp('===========================================')
disp('                DRIVER HYP                 ')
disp('===========================================')
disp('                                           ')
disp('     constitutive driver for numerical     ')
disp('        integration of hypoplastic         ')
disp('           constitutive equations          ')
disp('       using the algorithm proposed by     ')
disp('      Bardet & Choucair (1991), IJNAMG     ')
disp('                                           ')
disp('           written by C. Tamagnini         ')
disp('      modifications for hypoplasticity     ')
disp('                by D. Masin                ')
disp('                                           ')
disp('===========================================')
disp('                                           ')
disp('     ALERT Olek Zienkiewicz School 2012    ')
disp('      "Constitutive Modeling of Soils"     ')
disp('             Dresden, sept. 2012           ')
disp('                                           ')
disp('===========================================')

%% input data and initial material state

input_data;
init_state;

%% call integration procedure

[SS,EE,INV_S,INV_E,HARD] = update_1(y0,parms,nspb,path_info);

%% output section: writes results to various files and plot diagrams

save inv_s.txt INV_S -ascii;
save inv_e.txt INV_E -ascii;
save stresses.txt SS -ascii;
save strains.txt EE -ascii;
save hard.txt HARD -ascii;

xyplot_vmh;

