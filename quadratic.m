% Test bisection
clear
solver_flag = 1;
fun = @test_function02;

filter_list = [1e-15, 1e-2,1e-14,1e-2,2];
x_guess0 = 1;
guess_list1 = linspace(-10,-5,1000);
guess_list2 = linspace(2,7,1000);
    
convergence_analysis(solver_flag, fun, ...
x_guess0, guess_list1, guess_list2, filter_list)

hold off

function [f_val,dfdx] = test_function02(x)
    global input_list;
    input_list(:,end+1) = x;
    f_val = (x-30.879).^2;
    dfdx = 2*(x-30.879);
end
