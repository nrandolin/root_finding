function [f_val, dfdx] = test_function02(x)
    %declare input_list as a global variable
    global input_list;

    %append the current input to input_list
    %formatted so this works even if x is a column vector instead of a scalar
    input_list(:,end+1) = x;

    %perform the rest of the computation to generate output
    %I just put in a quadratic function as an example
    f_val = (x.^2) + 61.7580*(x) + 953.5126;
    dfdx = 2.*x + 61.7580;
end