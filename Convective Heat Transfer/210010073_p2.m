clear;
clc;
%2D Convection problem_ uniform slug flow in 1 direction
u = 1;
v = 0; %1D flow 
% Re Pr = uH/alpha  = 10*1 ; alpha =0.1;
H = 1;
L = 6;
imax =42;
jmax = 22;
dx= L/(imax-1);
dy = H/(jmax-1);
tc = dx/u;
a = 0.1;
dt = 0.1*min(0.5*dx^2*dy^2/(a*(dx^2 + dy^2)),tc);
%dt=0.001;

% dT/dt + udT/dx = a[d^2T/dx^6 + d^2T/dy^2]
thetas = zeros(jmax,imax);
thetas(:,1) = thetas(:,1) + 1;
%%{
thetas(:,imax) = 0*thetas(:,imax);
thetas(jmax,:) = 0*thetas(jmax,:);
thetas(1,:) = 0*thetas(1,:);
%}



%FOU SCHEME
e = 1000;
Tfo = thetas;
Tfn = Tfo;
while e >1e-4
    for i = 2:imax-1
        for j = 2:jmax-2
            E = a*dt/dx^2 - min(u,0)*dt/dx;
            W = a*dt/dx^2 + max(u,0)*dt/dx;
            N = a*dt/dy^2 - min(v,0)*dt/dy;
            S = a*dt/dy^2 + max(v,0)*dt/dy;
            b = (1-E-W-S-N)*Tfo(j,i);
            Tfn(j,i) = E*Tfo(j,i+1) + W*Tfo(j,i-1) + N*Tfo(j+1,i) +S*Tfo(j-1,i)+b;
        end
    end
    e = L^2*max(max(abs(Tfn - Tfo)))/(a*dt);
    Tfo = Tfn;
end

%CD SCHEME:
e =1;
Tco = thetas;
Tcn = thetas;
while e > 1e-4 
%while e < 2
    for i = 2:imax-1
        for j = 2:jmax-2
            Diffy = (Tco(j+1,i)+Tco(j-1,i) -2*Tco(j,i))/dy^2;
            Diffx = (Tco(j,i+1)+Tco(j,i-1) -2*Tco(j,i))/dx^2;
            adv  = u*(Tco(j,i+1) -Tco(j,i-1))/(2*dx) + v*(Tco(j+1,i) -Tco(j-1,i))/(2*dy);
            Tcn(j,i) = Tco(j,i) + dt*(a*Diffx + a*Diffy - adv);
        end
    end
    e = L^2*max(max(abs(Tcn - Tco)))/(a*dt);
    %e=e+1;
    Tco = Tcn;
end
%Symmetry enforcer
%{
for j = 1:jmax
    Tcn(j,:) =0.5*(Tcn(jmax+1-j,:) + Tcn(j,:));
    Tcn(jmax+1-j,:) = Tcn(j,:);
    Tfn(j,:) =0.5*(Tfn(jmax+1-j,:) + Tfn(j,:));
    Tfn(jmax+1-j,:) = Tfn(j,:);
end
%}


[x,y] =meshgrid(dx*(0:imax-1),dy*(0:jmax-1));
%%{
contourf(x,y,Tcn,'ShowText','on');
xlabel("x dist in m")
ylabel('y dist in m')
title("\theta distribution C.D")
colormap("parula");
cb = colorbar('eastoutside');
%}
%{
plot(Tcn(:,8),y(:,1),Tcn(:,15),y(:,1),Tcn(:,21),y(:,1),Tcn(:,35),y(:,1),'LineWidth',1.5)
legend('at x/H = 1','at x/H = 2','at x/H =3','at x/H =5');
ylabel("y distance in meters")
xlabel("\theta : non dimensional temperature")
title("\theta Profile along Y direction","Explicit C.D scheme at various x values")
grid on
%}
