clear;
clc;
%CHT Code for 2D plate conduction EXPLICIT METHOD
imax= 11;
%initial conditions:
T0 = 30;
T =zeros(imax,imax);
for i =1:imax
    for j=1:imax
        T(i,j)=T0;
    end
end
%Boundary conditions
Twb = 100;
Teb = 300;
Tsb = 200;
Tnb = 400;
for i =1:imax
    T(i,imax) = Teb;
    T(imax,i) = Tnb;
    T(i,1) = Twb;
    T(1,i) = Tsb;
end
%Dimentional stuff
L1 = 1;
L2 = 1;
dx = L1/(imax-1);
dy =L2/(imax-1);
Lc = L1;
est  = 1e-4;

%heat supply
Q  =100e3;
%Q=0;

%Material properties
Cp = 500;
rho = 7750;
k = 16.2;
dTc = Tnb -Twb;
a = k/(Cp*rho);
b = Q/(Cp*rho);
dt =0.1*dx*dx*dy*dy/(2*a*(dx*dx + dy*dy));
Tnew = T+0;
Told = T+0;
e = 10000;
%js=1;
while e > est
    for i = 2:imax-1
        for j = 2:imax-1
            D2Tx = (Told(i+1,j) + Told(i-1,j) -2*Told(i,j))/(dx*dx);
            D2Ty = (Told(i,j+1) + Told(i,j-1) -2*Told(i,j))/(dy*dy);
            Tnew(i,j) = Told(i,j) + a*dt*(D2Ty + D2Tx) + b*dt;
        end
    end
    for i =1:imax
    Tnew(i,imax) = Teb;
    Tnew(imax,i) = Tnb;
    Tnew(i,1) = Twb;
    Tnew(1,i) = Tsb;
    end
    e =(Lc*Lc/(a*dTc))*max(max(abs((Tnew-Told)/dt)));
    Told = Tnew;
    %Tp(js) = Tnew(5,5);
    %js = js+1;
end
Tnew(imax,imax) =0.5*(Tnew(imax,imax-1)+Tnew(imax-1,imax));
Tnew(1,imax) =0.5*(Tnew(1,imax-1)+Tnew(2,imax));
Tnew(imax,1) =0.5*(Tnew(imax,2)+Tnew(imax-1,1));
Tnew(1,1) =0.5*(Tnew(1,2)+Tnew(2,1));

[X,Y]= meshgrid(dx*(0:imax-1),dy*(0:imax-1));
contourf(X,Y,Tnew,'ShowText','on');
cb = colorbar('eastoutside');
colormap("turbo");
title("Temp distribution With heat generation")
