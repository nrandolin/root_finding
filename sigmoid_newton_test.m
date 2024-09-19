clear all

% test newton
%define sigmoid
a = 27.3; b = 2; c = 8.3; d = -3;
test_pair = @(x)[c*(exp((x-a)/b))/(1+ exp((x-a)/b))+d, c*((1+ exp((x-a)/b))*((exp((x-a)/b))/b)-(exp((x-a)/b))*((exp((x-a)/b))/b))./((1+ exp((x-a)/b))^2)];


test_newt = newton_brr(test_pair, 30.1);

low_x = linspace(0,22.2,1000);
mid_x = linspace(22.2,30.15,1000);
end_x = linspace(30.15,50, 1000);
low_y = test_function03(low_x);
mid_y = test_function03(mid_x);
end_y = test_function03(end_x);
line_x = linspace(0,50,100);
line_y = zeros(1,100);

figure()
plot(low_x, low_y, "color", "r")
axis([0,50, -8, 8])
hold on
plot(mid_x, mid_y, "Color", "b")
plot(end_x, end_y, "Color", "r")
plot(line_x, line_y, "Color", "black")
title("Newton Limits")
xlabel("x")
ylabel("f(x)")
hold off

clear
%% Test Secant and fsolve
a = 27.3; b = 2; c = 8.3; d = -3;
func = @(x) c*(exp((x-a)/b))/(1+ exp((x-a)/b))+d;
upper_limit = 32.2
root = secant_brr(func,upper_limit,32.3)
lower_limit = root - abs(root - upper_limit +0.1)

low_x = linspace(0,lower_limit,1000);
mid_x = linspace(lower_limit,upper_limit,1000);
end_x = linspace(upper_limit,50, 1000);
low_y = test_function03(low_x);
mid_y = test_function03(mid_x);
end_y = test_function03(end_x);
line_x = linspace(0,50,100);
line_y = zeros(1,100);

figure()
plot(low_x, low_y, "color", "r")
axis([0,50, -8, 8])
hold on
plot(mid_x, mid_y, "Color", "b")
plot(end_x, end_y, "Color", "r")
plot(line_x, line_y, "Color", "black")
title("Secant Limits")
xlabel("x")
ylabel("f(x)")
hold off

%% Test Fsolve
clear all
a = 27.3; b = 2; c = 8.3; d = -3;
func = @(x) c*(exp((x-a)/b))/(1+ exp((x-a)/b))+d;
upper_limit = 26
root = fsolve(func,upper_limit)
lower_limit = root - abs(root - upper_limit +0.1)

low_x = linspace(0,50,1000);
low_y = test_function03(low_x);
line_x = linspace(0,50,100);
line_y = zeros(1,100);

figure()
plot(low_x, low_y, "color", "r")
axis([0,50, -8, 8])
hold on
plot(line_x, line_y, "Color", "black")
title("Fsolve Limits (It has none)")
xlabel("x")
ylabel("f(x)")
hold off
%%
%Note that fun(x) should output [f,dfdx], where dfdx is the derivative of f
function z = newton_brr(fun,x0)
    difference = 1;
    x1 = x0;
    while difference > 10^-14
        x0 = x1;
        y_val = fun(x0); 
        if abs(y_val(2)) < 10^-14
            quit
        else
            x1 = x0 - y_val(1)/ y_val(2);
            difference = abs(x1 - x0);
        end 
    end
    z = x1;
end

function x = secant_brr(fun, x0, x1)
    y0 = fun(x0);
    y1 = fun(x1);
    difference = 1;
    while difference > 10^-14
        if abs(y1-y0) < 10^-14
            quit
        else
        x2 = x1-y1*((x1-x0)/(y1-y0));
        difference = abs(x2-x1);
        
        x0 = x1;
        y0 = y1;

        x1 = x2;
        y1 = fun(x1);
        end
    end
    x = x2;
end