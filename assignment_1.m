%Definition of the test function and its derivative
test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
test_derivative01 = @(x) 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;

test_func_1 = @(x)[(x.^2) - 9, (2*x)];
test_func = @(x) (x.^2) - 9;

zeros = bisection_solver(test_func, 0, 6);
test = bisection_solver(test_func01, 0, 10);
test_1 = newton_solver(test_func_1, 4)

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
        y_val = fun(x0) ;
        x1 = x0 - y_val(1)/ y_val(2);
        difference = abs(x1 - x0);
    end
    z = x1;
end



function x = secant_solver(fun, x0, x1)
    difference = 1;
    while difference > 10^-14
        x2 = x1-fun(x1)*((x1-x0)/(fun(x1)-fun(x0)));
        difference = abs(x2-x1);
        x0 = x1;
    end
    x = x2;
end

