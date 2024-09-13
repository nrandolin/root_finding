function [f_val,dfdx] = test_function02(x)
    global input_list;
    input_list(:,end+1) = x;
    f_val = (x-30.879).^2;
    dfdx = 2*(x-30.879);
end