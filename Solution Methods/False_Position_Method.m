%False Position method
clear all;
%intial guess
a = -1;
b =  0;
if abs(f(a))>abs(f(b))
    cold = b;
else
    cold =a;
end
c = (f(a)*b - f(b)*a)/(f(a) -f(b));
%Tolerance
e = 0.001;
%Loop
while (abs(c-cold) > e)
    if f(b)*f(c)<=0 
        a = c + 0.0;
    else 
        b = c + 0.0;
    end
    cold = c +0 ;
    c = (f(a)*b - f(b)*a)/(f(a) -f(b));
end
disp(c);


function p = f(x)
p = x^2 -2*x - 3;
end