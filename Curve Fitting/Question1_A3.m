%clear all;
clc;
%Code to solve the given specific problem
M = [ 2 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 2];
Ts = [0;0;0;0]; %[T2,T3,T4,T5]
Null = Ts + 0;
T1 = 0;
T6 = 100;
Tb = [ T1;0;0;T6];
%TDMA Solution:
Me = M + 0; 
Tbe = Tb + 0;% Lower part elimination code
for i = 2:4
    k = Me(i,i-1)/Me(i-1,i-1);
    Tbe(i) = Tbe(i) - k*Tbe(i-1);
    for j =1:4
        Me(i,j) = Me(i,j) - k*Me(i-1,j);
    end
end
%Back substitution and solving
Ts(4) = Tbe(4)/Me(4,4);
for i = 1:3
    Ts(4-i) = (Tbe(4-i) - Ts(5-i)*Me(4-i,5-i))/Me(4-i,4-i);
end
disp("TDMA Solution");
disp(Ts);
% ____ITERATIVE METHODS_______
% Initial Guesses
Ts0 = Null + 0;
Ts50 = Null + 50;
Ts100 = Null + 100;
% ----Jacobi Method ----
e = 100;
old  = Ts0;
n0 = 0;
n50 = 0;
n100 =0;
while e >= 1e-4
    n0 = n0 +1;
    new  = Null + 0;
    s = 0;
    for i = 1:4
        for j = 1:4
             s = s - M(i,j)*old(j);
        end
        new(i) = (Tb(i) + s + M(i,i)*old(i))/M(i,i);
    end
    e = max(abs(new - old));
    old = new + 0;
end
disp("Jacobi method with Initial iterates  (0,0,0,0)");
disp(new);
%------------------------------------------------
e = 100;
old  = Ts50;
while e >= 1e-4
    n50 = n50 + 1;
    new  = Null + 0;
    s = 0;
    for i = 1:4
        for j = 1:4
             s = s - M(i,j)*old(j);
        end
        new(i) = (Tb(i) + s + M(i,i)*old(i))/M(i,i);
    end
    e = max(abs(new - old));
    old = new + 0;
end
disp("Jacobi method with Initial iterates  (50,50,50,50)");
disp(new);

%-----------------------------------------
e = 100;
old  = Ts100;
while e >= 1e-4
    n100 = n100 +1;
    new  = Null + 0;
    s = 0;
    for i = 1:4
        for j = 1:4
             s = s - M(i,j)*old(j);
        end
        new(i) = (Tb(i) + s + M(i,i)*old(i))/M(i,i);
    end
    e = max(abs(new - old));
    old = new + 0;
end
disp("Jacobi method with Initial iterates  (100,100,100,100)");
disp(new);

%Comparative performance_______________________________________
disp("Number of iterations for 0  50 and 100 in Jacobi Method");
disp([n0 n50 n100]);
%GAUSS SIDEL METHOD______________________________________________
%Intial condition 50
e = 100;
old  = Ts50;
new = old +0;
n_50 =0;
while e > 1e-4
    n_50 = n_50 +1;
    for i = 1:4
        sb = 0;
        sa = 0;
        for j = 1:4
            if j < i
                sb = sb - new(j)*M(i,j);
            end
            if j > i
                sa = sa - old(j)*M(i,j);
            end
        end
        new(i) = (Tb(i) + sa +sb)/M(i,i);
    end
        
    e = max(abs(new - old));
    old = new;

end
disp("Gauss Sidel method with Initial iterates  (50,50,50,50)");
disp(new);

%Intial condition 0
e = 100;
old  = Ts0;
new = old +0;
n_0 =0;
while e > 1e-4
    n_0 = n_0 +1;
    for i = 1:4
        sb = 0;
        sa = 0;
        for j = 1:4
            if j < i
                sb = sb - new(j)*M(i,j);
            end
            if j > i
                sa = sa - old(j)*M(i,j);
            end
        end
        new(i) = (Tb(i) + sa +sb)/M(i,i);
    end
        
    e = max(abs(new - old));
    old = new;

end
disp("Gauss Sidel method with Initial iterates  (0,0,0,0)");
disp(new);

%Intial condition 100
e = 100;
old  = Ts100;
new = old +0;
n_100 =0;
while e > 1e-4
    n_100 = n_100 +1;
    for i = 1:4
        sb = 0;
        sa = 0;
        for j = 1:4
            if j < i
                sb = sb - new(j)*M(i,j);
            end
            if j > i
                sa = sa - old(j)*M(i,j);
            end
        end
        new(i) = (Tb(i) + sa +sb)/M(i,i);
    end
        
    e = max(abs(new - old));
    old = new;

end
disp("Gauss Sidel method with Initial iterates  (100,100,100,100)");
disp(new);

%Comparative performance_______________________________________
disp("Number of iterations for 0  50 and 100 in GS Method");
disp([n_0 n_50 n_100]);
