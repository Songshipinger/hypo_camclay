% plotting section
%
% sand hypoplastic constitutive model

array_sizes = size(INV_E);
nrow = array_sizes(1);

figure(1)
clf

% epss-q plot

subplot(2,2,1)

xx=INV_E(1:nrow,2);
yy=INV_S(1:nrow,2);

plot(xx,yy,'r-')
xlabel('deviatoric strain')
ylabel('deviator stress q [kPa]')
grid on

% p:q plot

subplot(2,2,2)

xx=INV_S(1:nrow,1);
yy=INV_S(1:nrow,2);

plot(xx,yy,'r-')
xlabel('mean effective stress p [kPa]')
ylabel('deviator stress q [kPa]')
grid on

% epsv-epss plot

subplot(2,2,3)

xx=INV_E(1:nrow,2);
yy=INV_E(1:nrow,1);
plot(xx,yy,'r-')

ylabel('volumetric strain')
xlabel('deviatoric strain')
grid on

% eps_v-p plot

subplot(2,2,4)

xx=INV_S(1:nrow,1);
yy=INV_E(1:nrow,1);

plot(xx,yy,'r-')
xlabel('mean effective stress p [kPa]')
ylabel('volumetric strain')
grid on

figure(2)
clf

% c-eps_s plot

xx=INV_E(1:nrow,2);
yy=HARD(1:nrow,1);

plot(xx,yy,'r-')
ylabel('void ratio [-]')
xlabel('deviatoric strain')
grid on