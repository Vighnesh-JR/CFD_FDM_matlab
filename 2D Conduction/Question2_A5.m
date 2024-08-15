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
%material and environmental properties:
qw = 10e3; %east boundary
Tinf = 30;
Twb = 100;
L1 =1;
L2=1;
h=100;
k=16.2;
Cp = 500;
rho = 7750;

%Heat generation
%Q=0;
Q=50e3;


%iteration and grid vals:
dx = L1/(imax-1);
dy = L2/(imax-1);
est = 1e-4;
egs = 1e-4;
dTc = Twb -Tinf;
Lc = L1;
a = k/(Cp*rho);
b = Q/(Cp*rho);
dt =0.1*dx*dx*dy*dy/(2*a*(dx*dx + dy*dy));


Told = T;
er = 1000;

%BC application:
for i = 1:imax
    Told(i,1)=Twb;
    Told(1,i)=Told(2,i);
    Told(i,imax)= Told(i,imax-1) - qw*dx/k;
    Told(imax,i) = (k*Told(imax-1,i)/dx + h*Tinf)/(k/dx +h);
end
while er > est
    Tint = Told;
    eg = 1000;
    Tref = Tint;
    while eg > egs
        for i=2:imax-1
            for j=2:imax-1
                Dtx = a*dt*(Tint(i+1,j)+Tint(i-1,j))/dx^2;
                Dty = a*dt*(Tint(i,j+1)+Tint(i,j-1))/dy^2;
                f = 1 + 2*dt*a*(1/dx^2 + 1/dy^2);
                Tint(i,j) = (Told(i,j) + Dtx + Dty + b*dt)/f;
            end
        end
        for i = 1:imax
            Tint(i,1)=Twb;
            Tint(1,i)=Tint(2,i);
            Tint(i,imax)= Tint(i,imax-1) + qw*dx/k;
            Tint(imax,i) = (k*Tint(imax-1,i)/dy + h*Tinf)/(k/dy +h);
        end
        eg = max(max(abs(Tref -Tint)));
        Tref = Tint;
    end
 
    er = (Lc*Lc/(a*dTc))*max(max(abs(Told -Tint)/dt));
    Told = Tint;
end

[X,Y]= meshgrid(dx*(0:imax-1),dy*(0:imax-1));
contourf(X,Y,Told,'ShowText','on');
colormap(summer);
cb = colorbar('eastoutside');
title("Temp distribution With heat generation")

