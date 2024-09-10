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

