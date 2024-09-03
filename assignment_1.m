%Definition of the test function and its derivative
test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
test_derivative01 = @(x) 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;
test_pair01 = @(x)[((x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6)), (3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6)];

test_func_1 = @(x)[(x.^2) - 9, (2*x)];
test_func = @(x) (x.^2) - 9;

test_bi = bisection_solver(test_func01, 0, 10);
% plot bisection
figure()
fplot(test_func01, [-10 35])
hold on
yline(0)
grid on
hold on
plot(test_bi, test_func01(test_bi), 'd')
title("Bisection Root Solver Output")
hold off

test_newt = newton_solver(test_pair01, 0.717);
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

test_sec = secant_solver(test_func01, 1, 4);
% plot secant
figure()
fplot(test_func01, [-10 35])
hold on
yline(0)
grid on
hold on
plot(test_sec, test_func01(test_sec), 'd')
title("Secant Root Solver Output")
hold off


function x = bisection_solver(fun,x_left,x_right)
    x_mid = 10;
    c = 11;
    while abs(x_mid - c) > 10^-14
        if fun(x_left) * fun(x_right) > 0 
            quit
        else
            x_mid = (x_left + x_right) / 2;
            if fun(x_left) * fun(x_mid) < 0
                x_right = x_mid;
                c = (x_left + x_right) / 2;
            elseif fun(x_right) * fun(x_mid) < 0
                x_left = x_mid;
                c = (x_right + x_left) / 2;
            elseif fun(x_mid) == 0
                c = x_mid;
                break
            end
        end
    end
    x = c;
end


%Note that fun(x) should output [f,dfdx], where dfdx is the derivative of f
function z = newton_solver(fun,x0)
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



function x = secant_solver(fun, x0, x1)
    difference = 1;
    while difference > 10^-14
        if abs(fun(x1)-fun(x0)) < 10^-14
            quit
        else
        x2 = x1-fun(x1)*((x1-x0)/(fun(x1)-fun(x0)));
        difference = abs(x2-x1);
        x0 = x1;
        x1 = x2;
        end
    end
    x = x2;
end

