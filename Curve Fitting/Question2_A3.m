clc;
%CURVE FITTING_______________________________________________
Re = [50 60 70 80 90 100 110 120];
Cd = [4.8504 3.3168 2.7111 2.3742 2.1465 1.9854 1.8618 1.7649];
Nu = [1.2095 1.5727 1.8417 2.0619 2.2510 2.4185 2.5689 2.7072];

%Logarithmic dataset_________________________________________
Lg_Cd = log(Cd);
Lg_Nu = log(Nu);
%Finding values for Functional Fit:___________________________
%Cd
disp("Linear Quadratic and exponential for Cd vs Re")
LC= lin_reg(Re,Cd);
QC = quad_reg(Re,Cd);
EC = lin_reg(Re,Lg_Cd);

disp(LC);
disp(QC);
disp([exp(EC(1)) EC(2)]);
%Nu
disp("Linear Quadratic and exponential for Nu vs Re")
LN= lin_reg(Re,Nu);
QN = quad_reg(Re,Nu);
EN= lin_reg(Re,Lg_Nu);

disp(LN);
disp(QN);
disp([exp(EN(1)),EN(2)]);



disp("plots")
%{
plot(Re,Cd,"r^",'LineWidth',3);
hold on 
plot(Re,Re*LC(2)+LC(1),"b-",'LineWidth',1.5)
plot(Re,QC(3)*Re.^2+QC(2)*Re+QC(1),"k-",'LineWidth',1.5)
plot(Re,exp(EC(2)*Re + EC(1)),"g-",'LineWidth',1.5)
title("Cd vs Re plot")
xlabel('Reynolds Number') 
ylabel('Coefficient of Drag') 
legend('Actual','Linear Regression','Quadratic regression','exponential');
%}

%{
plot(Re,Nu,"r^",'LineWidth',3);
hold on 
plot(Re,Re*LN(2)+LN(1),"b-",'LineWidth',1.5)
plot(Re,QN(3)*Re.^2+QN(2)*Re+QN(1),"k-",'LineWidth',1.5)
plot(Re,exp(EN(2)*Re + EN(1)),"g-",'LineWidth',1.5)
title("Nu vs Re plot")
xlabel('Reynolds Number') 
ylabel('Nusselt number ') 
legend('Actual','Linear Regression','Quadratic regression','exponential','Location','north');
%}




%Linear Regression algorithm:_________________________________
function L = lin_reg(x,y)
Sxx = mean(x.*x);
Sxy = mean(x.*y);
Xb = mean(x);
Yb = mean(y);
L(2) = (Sxy - Xb*Yb)/(Sxx - Xb*Xb);
L(1) = Yb - L(2)*Xb;
end
%Quadratic Regression algorithm
function Q =quad_reg(x,y)
A = [1,mean(x),mean(x.^2);mean(x),mean(x.^2),mean(x.^3);mean(x.^2),mean(x.^3),mean(x.^4)];
b = [mean(y);mean(y.*x);mean(y.*x.*x)];
Q = A\b;
Q = Q.';
end
