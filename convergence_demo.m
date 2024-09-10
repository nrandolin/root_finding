% Test bisection
clear
solver_flag = 1;
fun = @test_function;

filter_list = [1e-15, 1e-2,1e-14,1e-2,2];
x_guess0 = 1;
guess_list1 = linspace(-10,-5,1000);
guess_list2 = linspace(2,7,1000);
    
convergence_analysis(solver_flag, fun, ...
x_guess0, guess_list1, guess_list2, filter_list)

hold off



