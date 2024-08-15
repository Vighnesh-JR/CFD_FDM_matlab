clear;
clc;
L =0.01;
K= 16.2;
Q =100e6;
%Q=0;
cp=500;
p=7750;
imax = 21;
Twb =0;
a = K/(cp*p);
b = Q/(cp*p);
dx = L/(imax-1);
Xs=dx*(0:imax-1); 
dt = 0.1*(0.5*dx*dx/a);
h = 1000;
T0 = 100;
%Initial condition_____________________________
Ti = 30;
Ts = zeros(1,imax)+ Ti;
Ts(1) = Twb;
Told =Ts;
Tnew =Ts;
T_t = Ts;
Told(imax) = (h*T0 + K*Told(imax-1)/dx)/(h + K/dx);
e=1000;
while e>1e-3
    Tnew(1)=Twb;
    for i = 2:imax-1
        Tnew(i)=Told(i) + dt*(a*(Told(i+1)+Told(i-1)-2*Told(i))/(dx*dx) + b);
    end
    Tnew(imax) = (h*T0 + K*Tnew(imax-1)/dx)/(h + K/dx);
    e = max(abs(Tnew -Told));
    %Td = Told;
    Told = Tnew;
end
T_ref = -0.5*b*(Xs.^2)/a + Xs*((K*b*L/a) + 0.5*b*L*L*h/a + h*T0+h*Twb)/(K+h*L) +Twb;
%plotting code
plot(Xs,Tnew,"ro-",'Linewidth',2)
hold on
plot(Xs,T_ref,"k-",'Linewidth',1)
%title("Temperature distribution without heat generation")
xlabel("Thickness in m")
ylabel("Temperature in Celcius");
legend("Numerical","Analytical",'Location','south');