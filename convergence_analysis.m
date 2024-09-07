%% Function Start
%Example template for analysis function
%INPUTS:
%solver_flag: an integer from 1-4 indicating which solver to use
% 1->Bisection 2-> Newton 3->Secant 4->fzero
%fun: the mathematical function that we are using the
% solver to compute the root of
%x_guess0: the initial guess used to compute x_root
%guess_list1: a list of initial guesses for each trial
%guess_list2: a second list of initial guesses for each trial
% if guess_list2 is not needed, then set to zero in input
%filter_list: a list of constants used to filter the collected data
function convergence_analysis(solver_flag, fun, ...
x_guess0, guess_list1, guess_list2, filter_list)
    % Find actual* Root
    x_root = fzero(fun, x_guess0);

    % Check method and calculate actual p & K values
    if solver_flag == 1
        p_actual = 1;
    elseif solver_flag == 2
        p_actual = 2;
        [dfdx,d2fdx2] = approximate_derivative(fun,x_root);
        k_actual = d2fdx2/(2*dfdx);
    elseif solver_flag == 3
        p_actual = (1 +sqrt(5))/2;
        k_actual = (1 +sqrt(5))/2;
    elseif solver_flag == 4
        p_actual = "fzero";
        k_actual = "fzero";
    else
        print "Unknown method"
        return
    end

    
end

%% Derivative Function
function [dfdx,d2fdx2] = approximate_derivative(fun,x)
    %set the step size to be tiny
    delta_x = 1e-6;
    %compute the function at different points near x
    f_left = fun(x-delta_x);
    f_0 = fun(x);
    f_right = fun(x+delta_x);
    %approximate the first derivative
    dfdx = (f_right-f_left)/(2*delta_x);
    %approximate the second derivative
    d2fdx2 = (f_right-2*f_0+f_left)/(delta_x^2);
end
%% Linear Regression Function
%example for how to compute the fit line
%data points to be used in the regression
%x_regression -> e_n
%y_regression -> e_{n+1}
%p and k are the output coefficients
function [p,k] = generate_error_fit(x_regression,y_regression)
    %generate Y, X1, and X2
    %note that I use the transpose operator (’)
    %to convert the result from a row vector to a column
    %If you are copy-pasting, the ’ character may not work correctly
    Y = log(y_regression)';
    X1 = log(x_regression)';
    X2 = ones(length(X1),1);
    %run the regression
    coeff_vec = regress(Y,[X1,X2]);
    %pull out the coefficients from the fit
    p = coeff_vec(1);
    k = exp(coeff_vec(2));
end
%% Collect Data
function output = test_function(x)
    %declare input_list as a global variable
    global input_list;

    %append the current input to input_list
    %formatted so this works even if x is a column vector instead of a scalar
    input_list(:,end+1) = x;

    %perform the rest of the computation to generate output
    %I just put in a quadratic function as an example
    output = (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
end
%% bisection solver
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
%% newton solver
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


%% secant solver
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