clear all;
clc;
thetas = load("Thetas.txt");
%thetas = [1 1 1 1 1 1 1 1 1; 
          %0.489 0.6749 0.7546 0.7892 0.7990 0.7891 0.7546 0.6749 0.489;
          %0.2811 0.4562 0.5541 0.603 0.6179 0.603 0.554 0.4561 0.281;
          %0.1791 0.3145 0.4026 0.451 0.4663 0.4509 0.4025 0.3144 0.179;
          %0.1207 0.2202 0.2907 0.3318 0.3452 0.3317 0.2906 0.2201 0.1205];




%Numerical Differentiation
tw = 10; %upper wall temperature
t0 = 30; %Temperature of all other walls which is room temperature
Ts = thetas*(tw - t0);
k = 237; %Conductivity
deltax =0.1;
deltay = 0.1;
Flux_1 = zeros(1,9); %I-Order derivative
Flux_2 = zeros(1,9); %II-Order derivative
Flux_3 = zeros(1,9); %III-Order derivative
for i = 1:9
    Flux_1(i) =- k*(Ts(1,i)-Ts(2,i))/deltay;
    Flux_2(i) =- k*(3*Ts(1,i) - 4*Ts(2,i)+ Ts(3,i))/(2*deltay);
    Flux_3(i) =- k*(11*Ts(1,i)-18*Ts(2,i)+9*Ts(3,i)-2*Ts(4,i))/(6*deltay);
end
distances = deltax*(1:9);
plot(distances,Flux_1,'bs-',distances,Flux_2,'go-',distances,Flux_3,'rd-','LineWidth',2);
title("Heat Flux vs Distance on the wall");
xlabel("distance-x");
ylabel("Heat Flux");
legend("First order","Second Order","Third Order");



%Numerical Integration
Heat_flux = Flux_2 + 0;
Total_Heat_flow_1 = 0; %Trapezoidal Rule
Total_Heat_flow_2 = 0; %Simpsons 1/3rd Rule
for i = 1:8
    Total_Heat_flow_1 = Total_Heat_flow_1 + deltax*(Heat_flux(i)+Heat_flux(i+1))/2;
end
for i = 1:7
    Total_Heat_flow_2 = Total_Heat_flow_2 + deltax*(Heat_flux(i)+4*Heat_flux(i+1)+Heat_flux(i+2))/3;
end
disp(["Trapezoidal = ",Total_Heat_flow_1; " Simpsons 1/3rd rule ",Total_Heat_flow_2]);



