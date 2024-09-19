## Applied Math: Assignment 01 - Root Finding

An assignment exploring root-finding methods, including a problem with predicting the collision of an egg when on a given trajectory.

## Root Finding
To run the convergence analysis functions with the original function, open and run convergence_demo.m, you can comment in/out the different sections for each solver type. 
To run the convergence analysis functions with a quadratic we made up, open and run convergence_test01.m, you can comment in/out the different sections for each solver type. 
To run the convergence analysis functions with the quadratic whoos root is at the minimumn, open and run convergence_test02.m, you can comment in/out the different sections for each solver type. 
To run the convergence analysis functions with sigmoid function, open and run convergence_test03, you can comment in/out the different sections for each solver type. 
## Egg Collision
To run the collision prediction problem for a tumbling egg, run *egg.m*. Upon running this file, an animation of the egg collision plotting will be output, which will be saved as a video. If the user runs this file, they will need to change the file path noted in the function 
*animation* to a file path on their computer. 
## File guide
approximate_derivative.m approximates the derivative
assignment_1.m contains the original code for creating the root finding algorithums
bisection_solver.m is a function for the bisection solver
convergence_analysis.m is a function that creates the convergence analysis plots, and calculates p/k for each function.
generate_error_fit.m is a function that generates a line of best fit
global_bisection.m is a function that runs the bisection solver with input_list as a global variable
global_newton.m is a function that runs the newton solver with input_list as a global variable
global_secant.m is a function thatruns the secant solver with input_list as a global variable
newton_solver.m is a function thatfor the secant
secant_solver.m is a function for the secant solver
sigmoid_newton_test.m tests the newton and secant functions with the sigmoid graph to find the convergence interval
test_function.m is a function containing the default equasion/function
test_function01.m is a function containing an arbitrary test equ/funtion
test_function02.m is a function containing the quadratic with a root at the minimum 
test_function03.m is a function containing the sigmoid function
