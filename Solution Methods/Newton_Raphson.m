%Newton Raphson Method
clear all;
%intial condition
xold = -0.5;
%Iterate
x = xold - f(xold)/df(xold); 
% Tolerance
e =0.001;
while (abs(x-xold)>e)
    xold = x;
    x = xold - f(xold)/df(xold);
end
disp(x);

function p = f(x)
p = x^2 -2*x - 3;
end
function q = df(x)
q = 2*x -2;
end
