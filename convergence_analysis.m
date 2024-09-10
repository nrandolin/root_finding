clear
solver_flag = 1;
fun = @test_function;

x_guess0 =1;
guess_list1 = linspace(-10,-5,1000);
guess_list2 = linspace(5,25,1000);
filter_list = 1;

convergence_anal(solver_flag, fun, ...
x_guess0, guess_list1, guess_list2, filter_list)


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
function convergence_anal(solver_flag, fun, ...
x_guess0, guess_list1, guess_list2, filter_list)
    % Find actual* Root
    x_root = fzero(fun, x_guess0);

% bisection
    if solver_flag == 1
        % Check method and calculate actual p & K values
        p_actual = 1;
        
        % do the global variable ing
        global input_list;
        %number of trials we would like to perform
        num_iter = length(guess_list1);
        %list for the left and right guesses that we would like
        %to use each trial. These guesses have all been chosen
        %so that each trial will converge to the same root
        %because the root is somewhere between -5 and 5.
        x_left_list = guess_list1;
        x_right_list = guess_list2;
        %list of estimate at current iteration (x_{n})
        %compiled across all trials
        x_current_list = [];
        %list of estimate at next iteration (x_{n+1})
        %compiled across all trials
        x_next_list = [];
        %keeps track of which iteration (n) in a trial
        %each data point was collected from
        index_list = [];
        %loop through each trial
        for n = 1:num_iter
            %pull out the left and right guess for the trial
            x_left = x_left_list(n);
            x_right = x_right_list(n);
            %clear the input_list global variable
            input_list = [];
            %run the bisection solver
            global_bisection(fun, x_left, x_right)
            %at this point, input_list will be populated with the values that
            %the solver called at each iteration.
            %In other words, it is now [x_1,x_2,...x_n-1,x_n]
            %append the collected data to the compilation
            x_current_list = [x_current_list,input_list(1:end-1)];
            x_next_list = [x_next_list,input_list(2:end)];
            index_list = [index_list,1:length(input_list)-1];
        end
        %At this point, x_current_list corresponds to many many
        %measurements of x_{n} across many trials
        %and x_next_list corresponds to many many measurements of
        %the corresponding value of x_{n+1} across many trials
        %this is the data the you want to clean and analaze

        % Newton
    elseif solver_flag == 2
        % actual p and k
        p_actual = 2;
        [dfdx, d2fdx2] = approximate_derivative(fun,x_root);
        k_actual = d2fdx2/(2*dfdx);
% Do the global variable thing
        num_iter = length(guess_list1);
        for n = 1:num_iter
            %pull out the left and right guess for the trial
            x_left = guess_list1(n);
            %clear the input_list global variable
            input_list = [];
            %run the bisection solver
            global_newton(test_function(fun), x_left)
            %at this point, input_list will be populated with the values that
            %the solver called at each iteration.
            %In other words, it is now [x_1,x_2,...x_n-1,x_n]
            %append the collected data to the compilation
            x_current_list = [x_current_list,input_list(1:end-1)];
            x_next_list = [x_next_list,input_list(2:end)];
            index_list = [index_list,1:length(input_list)-1];
        end



        % Secant
    elseif solver_flag == 3
        p_actual = (1 +sqrt(5))/2;
        k_actual = (1 +sqrt(5))/2;
        num_iter = length(guess_list1);

        % Do the global variable thing
        for n = 1:num_iter
            %pull out the left and right guess for the trial
            x_left = guess_list1(n);
            x_right = guess_list2(n);
            %clear the input_list global variable
            input_list = [];
            %run the bisection solver
            global_secant(test_function(fun), x_left, x_right)
            %at this point, input_list will be populated with the values that
            %the solver called at each iteration.
            %In other words, it is now [x_1,x_2,...x_n-1,x_n]
            %append the collected data to the compilation
            x_current_list = [x_current_list,input_list(1:end-1)];
            x_next_list = [x_next_list,input_list(2:end)];
            index_list = [index_list,1:length(input_list)-1];
        end
        % matlab
    elseif solver_flag == 4
        p_actual = "fzero";
        k_actual = "fzero";
    % Do the global variable thing
        num_iter = length(guess_list1);
        for n = 1:num_iter
            %pull out the left and right guess for the trial
            x_left = guess_list1(n);
            %clear the input_list global variable
            input_list = [];
            %run the bisection solver
            global_newton(test_function(fun), x_left)
            %at this point, input_list will be populated with the values that
            %the solver called at each iteration.
            %In other words, it is now [x_1,x_2,...x_n-1,x_n]
            %append the collected data to the compilation
            x_current_list = [x_current_list,input_list(1:end-1)];
            x_next_list = [x_next_list,input_list(2:end)];
            index_list = [index_list,1:length(input_list)-1];
        end
    else
        print "Unknown method"
        return
    end

    error_list0 = abs(x_current_list - x_root);
    error_list1 = abs(x_next_list - x_root);
    
    figure(1)
    loglog(error_list0, error_list1, 'ro', 'markersize', 1)
    
    %example for how to filter the error data
    %currently have error_list0, error_list1, index_list
    %data points to be used in the regression
    x_regression = []; % e_n
    y_regression = []; % e_{n+1}
    %iterate through the collected data
    for n=1:length(index_list)
        %if the error is not too big or too small
        %and it was enough iterations into the trial...
        if error_list0(n)>1e-15 && error_list0(n)<1e-2 && ...
        error_list1(n)>1e-14 && error_list1(n)<1e-2 && ...
        index_list(n)>2
            %then add it to the set of points for regression
            x_regression(end+1) = error_list0(n);
            y_regression(end+1) = error_list1(n);
        end
    end


    figure(2)
    loglog(x_regression, y_regression, 'ro', 'markersize', 1)
    hold on
    
    %example for how to plot fit line
    %generate x data on a logarithmic range
    [p,k] = generate_error_fit(x_regression,y_regression);
    fit_line_x = 10.^[-16:.01:1];
    %compute the corresponding y values
    fit_line_y = k*fit_line_x.^p;
    %plot on a loglog plot.
    loglog(fit_line_x,fit_line_y,'k-','linewidth',2)
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

%% Global Bisection
function global_bisection(func, x_left, x_right)
    %declare input_list as a global variable
    global input_list;

    %initialize input_list to be an empty array
    input_list = [];

    %initialize guesses for bisection solver

    %run the bisection solver
    x_root = bisection_solver(func,x_left,x_right);

    %at this point, input_list will be populated with the input arguments
    %that bisection_solver used to call test_function
    %plot the inputs
    % plot(1:length(input_list),input_list,'ko','markerfacecolor', 'k');

    %reset input_list for the next test
    %input_list = [];
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

%% Global Newton
function global_newton(func, x_0)
    %declare input_list as a global variable
    global input_list;

    %initialize input_list to be an empty array
    input_list = [];

    %initialize guesses for bisection solver

    %run the bisection solver
    x_root = newton_solver(func,x0);

    %at this point, input_list will be populated with the input arguments
    %that bisection_solver used to call test_function
    %plot the inputs
    % plot(1:length(input_list),input_list,'ko','markerfacecolor', 'k');

    %reset input_list for the next test
    %input_list = [];
end

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

%% Global Secant
function global_secant(fun, x0, x1)
    %declare input_list as a global variable
    global input_list;

    %initialize input_list to be an empty array
    input_list = [];

    %initialize guesses for bisection solver

    %run the bisection solver
    x_root = secant_solver(fun, x0, x1);

    %at this point, input_list will be populated with the input arguments
    %that bisection_solver used to call test_function
    %plot the inputs
    % plot(1:length(input_list),input_list,'ko','markerfacecolor', 'k');

    %reset input_list for the next test
    %input_list = [];
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
    %% Test Function
function [f_val, dfdx] = test_function(x)
    %declare input_list as a global variable
    global input_list;

    %append the current input to input_list
    %formatted so this works even if x is a column vector instead of a scalar
    input_list(:,end+1) = x;

    %perform the rest of the computation to generate output
    %I just put in a quadratic function as an example
    f_val = (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
    dfdx = 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;
end
