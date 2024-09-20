%% newton solver
function z = newton_solver(fun, x0)
    % set the difference to 1
       difference = 1;
    x1 = x0;
    %while the distance ist close enough to zero
    while difference > 10^-14
        x0 = x1;
        [y_val,d_val] = fun(x0);
        % make sure it dosent go to infinity
        if abs(d_val) < 10^-14
            quit
        % reset the point to the function of the tangent zero
        else
            x1 = x0 - y_val/ d_val;
            difference = abs(x1 - x0);
        end 
    end
    z = x1;
end
