% Secant Method
clear all;
% Initialisation
a = -0.5;
b = 0;
% Closeness to root
if abs(f(a)) > abs(f(b))
    xnm1 = b; % x(n -1)
    xnm2 = a; % x(n-2)
else
    xnm1 = a; % x(n -1)
    xnm2 = b; % x(n-2)
end
%The iterate
xn = (xnm1*f(xnm2) - f(xnm1)*xnm2)/(f(xnm2) - f(xnm1));
%tolerance
e =0.001;
while (abs(xn- xnm1)>e)
    xnm2 = xnm1 + 0;
    xnm1 = xn +0 ;
    xn = (xnm1*f(xnm2) - f(xnm1)*xnm2)/(f(xnm2) - f(xnm1));
end
disp(xn)

function p = f(x)
p = x^2 -2*x - 3;
end