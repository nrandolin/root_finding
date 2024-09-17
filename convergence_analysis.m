%% Function Start
function [p_predict, k_predict, p, k] = convergence_analysis(solver_flag, fun, ...
x_guess0, guess_list1, guess_list2, filter_list)
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

    % Define initial variables
    x_root = fzero(fun, x_guess0); % Known zero
    global input_list; % list of input values (error) from all iterations
    num_iter = length(guess_list1); % number of iterations
    x_left_list = guess_list1; % left guesses
    x_right_list = guess_list2; % right guesses (0 if newton)
    
    % Create list variables
    x_current_list = [];
    x_next_list = [];
    index_list = [];

    % Collect Error Data
    for n = 1:num_iter
        %pull out the left and right guess for the trial
        x_left = x_left_list(n);
        x_right = x_right_list(n);
        %clear the input_list global variable
        input_list = [];
        %run the solver
        if solver_flag == 1 % bisection
            p_predict = 1;
            k_predict = 1/2;
            global_bisection(fun, x_left, x_right)
        elseif solver_flag == 2 % newton
            p_predict = 2;
            [dfdx,d2fdx2] = approximate_derivative(fun,x_root);
            k_predict = d2fdx2/(2*dfdx);
            global_newton(fun, x_left)
        elseif solver_flag == 3 % secant
            p_predict = (1 +sqrt(5))/2;
            k_predict = (1 +sqrt(5))/2;
            global_secant(fun, x_left, x_left+0.1)
        elseif solver_flag == 4 % newton
            fzero(fun, x_left)
            p_predict = 0;
            k_predict = 0;
        else
            print "Unknown method"
            return
        end

        
        x_current_list = [x_current_list,input_list(1:end-1)];
        x_next_list = [x_next_list,input_list(2:end)];
        index_list = [index_list,1:length(input_list)-1];
    end

    % Calculte Error
    error_list0 = abs(x_current_list - x_root);
    error_list1 = abs(x_next_list - x_root);
    
    % Plot all Error data
    figure(1)
    loglog(error_list0, error_list1, 'ro', 'markersize', 1)
    title("Error Data")
    
    % Clean data
    x_regression = []; % e_n
    y_regression = []; % e_{n+1}
    %iterate through the collected data
    for n=1:length(index_list)
        %if the error is not too big or too small
        %and it was enough iterations into the trial...
        if error_list0(n)>filter_list(1) && error_list0(n)<filter_list(2) && ...
        error_list1(n)>filter_list(3) && error_list1(n)<filter_list(4) && ...
        index_list(n)>filter_list(5)
            %then add it to the set of points for regression
            x_regression(end+1) = error_list0(n);
            y_regression(end+1) = error_list1(n);
        end
    end

    % Plot clean data
    figure(2)
    loglog(x_regression, y_regression, 'ro', 'markersize', 1)
    hold on
    
    % pgenerate fit line
    [p,k] = generate_error_fit(x_regression,y_regression);
    fit_line_x = 10.^[-16:.01:1];
    %compute the corresponding y values
    fit_line_y = k*fit_line_x.^p;
    %plot on a loglog plot.
    loglog(fit_line_x,fit_line_y,'k-','linewidth',2)
end
