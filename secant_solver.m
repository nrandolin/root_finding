%% secant solver
function x = secant_solver(fun, x0, x1)
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