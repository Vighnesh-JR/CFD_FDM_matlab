% Bisection method
clear all;
%clc;
%intial guess
a = -1;
b =  0;
c = (a + b)/2;
%Tolerance
e = 0.001;
while abs(a-c) > e
    if f(a)*f(c)<=0 
        b = c + 0.0;
    else 
        a = c + 0.0;
    end
    c = (a+b)/2;
end
disp(c)


function p = f(x)
p = x^2 -2*x - 3;
end