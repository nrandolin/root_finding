function [f_val, dfdx] = test_function01(x)
    %declare input_list as a global variable
    global input_list;

    %append the current input to input_list
    %formatted so this works even if x is a column vector instead of a scalar
    input_list(:,end+1) = x;

    %perform the rest of the computation to generate output
    %I just put in a quadratic function as an example
    f_val = (x.^3) - 5.*(x.^2) + x - 3;
    dfdx = 3*(x.^2) - 10*x + 1;
end