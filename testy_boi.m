
fun = [test_function, test_function_deriv];
x0=1;
%Note that fun(x) should output [f,dfdx], where dfdx is the derivative of f
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

% plot newton
figure()
fplot(test_func01, [-10 35])
hold on
yline(0)
grid on
hold on
plot(test_newt, test_func01(test_newt), 'd')
title("Newton Root Solver Output")
hold off





