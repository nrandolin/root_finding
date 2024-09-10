%% newton solver
function z = newton_solver(fun, x0)
       difference = 1;
    x1 = x0;
    while difference > 10^-14
        x0 = x1;
        [y_val,d_val] = fun(x0); 
        if abs(d_val) < 10^-14
            quit
        else
            x1 = x0 - y_val/ d_val;
            difference = abs(x1 - x0);
        end 
    end
    z = x1;
end
