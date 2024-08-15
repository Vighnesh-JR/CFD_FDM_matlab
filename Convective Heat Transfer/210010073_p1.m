%PROBLEM 1 - 1D CONVECTION PROBLEM___________________
clear;
clc;

%Grid parameters:
imax = 11;
L=1;
dx = L/(imax-1);


%Material Properties
rho = 1;
u = 1;
Gamma = 0.02;

%Grid generation:
phis = zeros(imax,1);
Xs = dx*(0:imax-1);

e_max = 1e-4;
%boundary condition
phis(1)=0;
phis(imax)=1;

%FOU Scheme

Mf = zeros(imax-2,imax-2);
for i = 1:imax-2
    Mf(i,i) = 2*Gamma/dx^2 + rho*u/dx;
end
for i = 1:imax-3
    Mf(i,i+1) = rho*min(u,0)/dx - Gamma/dx^2;
    Mf(i+1,i) = rho*max(u,0)/dx - Gamma/dx^2;
end

b = zeros(imax-2,1);
b(1) = (-rho*max(u,0)/dx + Gamma/dx^2)*phis(1);
b(imax-2) = -(rho*min(u,0)/dx - Gamma/dx^2)*phis(imax);
phif = Mf\b;
phif = [phis(1);phif;phis(imax)];

%CD Scheme___
phic = phis;
%{
for i =2:imax-1
    phic(i) = 0.5*(phis(1) + phis(imax));
end
%}
phic_o = phic;
e=100;
a=0.1; %Relaxation factor
pe = rho*u*dx/Gamma;


while e>1e-3
    phic_o = phic;
    for i = 2:imax-1
        phic(i) = 0.5*((1-pe/2)*phic(i+1) + (1+pe/2)*phic(i-1)); 
    end
    e = max(abs(phic - phic_o));
    for i = 2:imax-1
        phic(i) = phic_o(i) + a*(phic(i) - phic_o(i));
    end
    

end
Phiref = phis(1) + (exp(Xs*rho*u/Gamma)-1)*(phis(imax) - phis(1))/(exp(rho*u/(Gamma*L))-1);


plot(Xs,phif,'o-',Xs,phic,'d-',Xs,Phiref,'k--','LineWidth',1.5);
legend('FOU scheme','CD scheme','Analytical','Location','northwest');
title('CONVECTED VARIABLE ACROSS L=1m','11 grid points')
xlabel('x');
ylabel('\phi')
grid on

