clear;
clc;
imax = 21;
L=0.01;
dx = L/(imax-1);
Mat = zeros(imax-2,imax-2);
for j = 1:imax-2
    Mat(j,j) = 2;
end
for j = 1:imax-3
    Mat(j,j+1)= -1;
end
for j = 1:imax-3
    Mat(j+1,j)=-1;
end
Twb = 0;
Teb = 100;
Tbnd = zeros(imax-2,1);
Tbnd(1)=Twb;
Tbnd(imax-2)=Teb;
Qgen = 0;
Qgen = 100e6; %MW/m^3;
k = 16.2; %W/m-K
B = Qgen/k;
Gen_vec = zeros(imax-2,1) + B*(dx^2);
Ddt = Gen_vec + Tbnd;
% Solution 
Ts = zeros(imax-2,1);
e = 100;
old  = Ts + 50; 
new = old +0;
while e > 1e-3
    for i = 1:imax-2
        sb = 0;
        sa = 0;
        for j = 1:imax-2
            if j < i
                sb = sb - new(j)*Mat(i,j);
            end
            if j > i
                sa = sa - old(j)*Mat(i,j);
            end
        end
        new(i) = (Ddt(i) + sa +sb)/Mat(i,i);
    end
        
    e = max(abs(new - old));
    old = new;
end
Ts = new;
Xs=dx*(0:imax-1);
Tc = [Twb;Ts;Teb];
T_ref = Twb + ((Teb - Twb)/L + B*L/2)*Xs - 0.5*B*(Xs.^2);
plot(Xs,Tc,"ro-",'Linewidth',2)
hold on
plot(Xs,T_ref,"k-",'Linewidth',1)
%title("Temperature distribution with heat generation")
xlabel("Thickness in m")
ylabel("Temperature in Celcius");
legend("Numerical","Analytical");
