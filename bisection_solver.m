%% bisection solver
function x = bisection_solver(fun,x_left,x_right)
    x_mid = 10;
    c = 11;
    % Define y compments for left and right
    y_left = fun(x_left);
    y_right = fun(x_right);
    while abs(x_mid - c) > 10^-14
        % If both are same sign, quit
        if y_left * y_right > 0 
            disp("your gess is wrong");
            x = 0;
            return 
        else
            % Define x mid and coresponding y
            x_mid = (x_left + x_right) / 2;
            y_mid = fun(x_mid);
            % If if yleft and y mid have different signs
            if y_left * y_mid < 0
                % move the midpoint to x_right
                x_right = x_mid;
                c = (x_left + x_right) / 2;
                % Recalculate y right to be at new right point
                y_right = y_mid;
            % If y right and y mid are different signs
            elseif y_right * y_mid < 0
                % move the midpoint to x_left
                x_left = x_mid;
                c = (x_right + x_left) / 2;
                % rcalculate y_left with new y point
                y_left = y_mid;
            elseif y_mid == 0
                c = x_mid;
                break
            end
        end
    end
    x = c;
end
