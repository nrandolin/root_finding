% % Test Newton
% clear
% solver_flag = 2;
% fun = @test_function02;
% 
% filter_list = [1e-15, 1e-2,1e-14,1e-2,2];
% x_guess0 = 1;
% guess_list1 = linspace(31,33,100);
% guess_list2 = linspace(2,7,100);
%     
% convergence_analysis(solver_flag, fun, ...
% x_guess0, guess_list1, guess_list2, filter_list)
% 
% hold off
%%
% % Test Secant
% clear
% solver_flag = 3;
% fun = @test_function02;
% 
% filter_list = [1e-15, 1e-2,1e-14,1e-2,2];
% x_guess0 = 1;
% guess_list1 = linspace(31,40,100);
% guess_list2 = linspace(40,49,100);
%     
% [p_precit, k_precidt, p, k] = convergence_analysis(solver_flag, fun, ...
% x_guess0, guess_list1, guess_list2, filter_list)
% 
% hold off
%%
% Test fzero
clear
solver_flag = 4;
fun = @test_function02;

filter_list = [1e-15, 1e-2,1e-14,1e-2,2];
x_guess0 = 1;
guess_list1 = linspace(31,32,100);
guess_list2 = linspace(40,49,100);
    
convergence_analysis(solver_flag, fun, ...
x_guess0, guess_list1, guess_list2, filter_list)

hold off


function [f_val, dfdx] = test_function02(x)
    %declare input_list as a global variable
    global input_list;

    %append the current input to input_list
    %formatted so this works even if x is a column vector instead of a scalar
    input_list(:,end+1) = x;

    %perform the rest of the computation to generate output
    %I just put in a quadratic function as an example
    f_val = (x-5).^2;
    dfdx = 2.*(x-5);
end

