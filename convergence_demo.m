% Test bisection
clear
solver_flag = 1;
fun = @test_function;

filter_list = [1e-15, 1e-2,1e-14,1e-2,2];
x_guess0 = 1;
guess_list1 = linspace(-10,-5,1000);
guess_list2 = linspace(2,7,1000);
    
[p_precit, k_precidt, p, k] = convergence_analysis(solver_flag, fun, ...
x_guess0, guess_list1, guess_list2, filter_list)

hold off

% clear
% 
% % Test Newton
% clear
% solver_flag = 2;
% fun = @test_function;
% 
% filter_list = [1e-15, 1e-2,1e-14,1e-2,2];
% x_guess0 = 1;
% guess_list1 = linspace(-1,3,100);
% guess_list2 = linspace(2,7,1000);
%     
% [p_precit, k_precidt, p, k] = convergence_analysis(solver_flag, fun, ...
% x_guess0, guess_list1, guess_list2, filter_list)
% 
% hold off

% % Test Secant
% clear
% solver_flag = 3;
% fun = @test_function;
% 
% filter_list = [1e-15, 1e-2,1e-14,1e-2,2];
% x_guess0 = 1;
% guess_list1 = linspace(-10,-5,1000);
% guess_list2 = linspace(2,7,1000);
%     
% [p_precit, k_precidt, p, k] = convergence_analysis(solver_flag, fun, ...
% x_guess0, guess_list1, guess_list2, filter_list)
% 
% hold off

% % Test fzero
% clear
% solver_flag = 4;
% fun = @test_function;
% 
% filter_list = [1e-15, 1e-2,1e-14,1e-2,2];
% x_guess0 = 1;
% guess_list1 = linspace(-2,3,1000);
% guess_list2 = linspace(2,7,1000);
%     
% [p_precit, k_precidt, p, k] = convergence_analysis(solver_flag, fun, ...
% x_guess0, guess_list1, guess_list2, filter_list)
% 
% hold off




