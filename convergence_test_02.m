% % Test Newton
% clear
% solver_flag = 2;
% fun = @test_function;
% 
% filter_list = [1e-15, 1e-2,1e-14,1e-2,2];
% x_guess0 = 1;
% guess_list1 = linspace(28,32,100);
% guess_list2 = linspace(2,7,100);
%     
% convergence_analysis(solver_flag, fun, ...
% x_guess0, guess_list1, guess_list2, filter_list)
% 
% hold off

% Test Secant
clear
solver_flag = 3;
fun = @test_function01;

filter_list = [1e-15, 1e-2,1e-14,1e-2,2];
x_guess0 = 1;
guess_list1 = linspace(31,40,100);
guess_list2 = linspace(40,49,100);
    
convergence_analysis(solver_flag, fun, ...
x_guess0, guess_list1, guess_list2, filter_list)

hold off

% % Test fzero
% clear
% solver_flag = 4;
% fun = @test_function01;
% 
% filter_list = [1e-15, 1e-2,1e-14,1e-2,2];
% x_guess0 = 1;
% guess_list1 = linspace(31,40,100);
% guess_list2 = linspace(40,49,100);
%     
% convergence_analysis(solver_flag, fun, ...
% x_guess0, guess_list1, guess_list2, filter_list)
% 
% hold off




