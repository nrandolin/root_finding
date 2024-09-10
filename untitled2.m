clear

solver_flag = 2;
fun = @test_function;

x_guess0 = 1;
guess_list1 = linspace(0,3,100);
guess_list2 = linspace(5,25,1000);
filter_list = 1;




function global_newton(func, x_0)
    %declare input_list as a global variable
    global input_list;

    %initialize input_list to be an empty array
    input_list = [];

    %initialize guesses for bisection solver

    %run the bisection solver
    x_root = newton_solver(func,x_0);

    %at this point, input_list will be populated with the input arguments
    %that bisection_solver used to call test_function
    %plot the inputs
    % plot(1:length(input_list),input_list,'ko','markerfacecolor', 'k');

    %reset input_list for the next test
    %input_list = [];
end


%% newton solver
%Note that fun(x) should output [f,dfdx], where dfdx is the derivative of f
function z = newton_solver(fun,x0)
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
