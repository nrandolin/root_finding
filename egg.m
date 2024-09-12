egg_params = struct();
egg_params.a = 3; egg_params.b = 2; egg_params.c = .15;
%specify the position and orientation of the egg
x0 = 5; y0 = 5; theta = pi/6;
[x_range,y_range] = compute_bounding_box(x0,y0,theta,egg_params);
figure;
eggxample01()
rectangle('Position',[x_range(1) y_range(1) x_range(2)-x_range(1) y_range(2)]);

%% WRAPPER 2
%set the oval hyper-parameters (MOVED TO BOUNDING BOX)
function [s_root_x, s_root_x_2, s_root_y, s_root_y_2] = wrapper_2()
    egg_params = struct();
    egg_params.a = 3; egg_params.b = 2; egg_params.c = .15;
    %specify the position and orientation of the egg
    x0 = 5; y0 = 5; theta = pi/6;
    %s = 0.7;
    %wrapper function that calls egg_wrapper1
    %but only takes s as an input (other inputs are fixed)
    %(single input)
    egg_wrapper_2x = @(s) egg_wrapper_x(s,x0,y0,theta,egg_params);
    %compute the value of s for which the corresponding point on the oval
    %has an x-coordinate of zero
    egg_wrapper_2y = @(s) egg_wrapper_y(s,x0,y0,theta,egg_params);
    s_root_x = secant_solver(egg_wrapper_2x,0,0.1);
    s_root_x_2 = secant_solver(egg_wrapper_2x,-1,0.1);
    s_root_y = secant_solver(egg_wrapper_2y,0,0.1);
    s_root_y_2 = secant_solver(egg_wrapper_2y,0,0.9);
end
%% EGG CALL
%template for how to properly call egg_func
%also provides example for how to interpret outputs
function eggxample01()
    %set the oval hyper-parameters
    egg_params = struct();
    egg_params.a = 3; egg_params.b = 2; egg_params.c = .15;
    %specify the position and orientation of the egg
    x0 = 5; y0 = 5; theta = pi/6;
    %set up the axis
    hold on; axis equal; axis square
    axis([0,10,0,10])
    %plot the origin of the egg frame
    plot(x0,y0,'ro','markerfacecolor','r');
    %compute the perimeter of the egg
    [V_list, G_list] = egg_func(linspace(0,1,100),x0,y0,theta,egg_params);
    %plot the perimeter of the egg
    plot(V_list(1,:),V_list(2,:),'k');
    %compute a single point along the egg (s=.8)
    %as well as the tangent vector at that point
    [V_single, G_single] = egg_func(.8,x0,y0,theta,egg_params);
    %plot this single point on the egg
    plot(V_single(1),V_single(2),'ro','markerfacecolor','r');
    %plot this tangent vector on the egg
    vector_scaling = .1;
    tan_vec_x = [V_single(1),V_single(1)+vector_scaling*G_single(1)];
    tan_vec_y = [V_single(2),V_single(2)+vector_scaling*G_single(2)];
    plot(tan_vec_x,tan_vec_y,'g')
    hold on;
end

%% BOUNDING BOX
%Function that computes the bounding box of an oval
%INPUTS:
%theta: rotation of the oval. theta is a number from 0 to 2*pi.
%x0: horizontal offset of the oval
%y0: vertical offset of the oval
%egg_params: a struct describing the hyperparameters of the oval
%OUTPUTS:
%x_range: the x limits of the bounding box in the form [x_min,x_max]
%y_range: the y limits of the bounding box in the form [y_min,y_max]
function [x_range,y_range] = compute_bounding_box(x0,y0,theta,egg_params)
    egg_wrapper_2x = @(s) egg_wrapper_x(s,x0,y0,theta,egg_params);
    %compute the value of s for which the corresponding point on the oval
    %has an x-coordinate of zero
    egg_wrapper_2y = @(s) egg_wrapper_y(s,x0,y0,theta,egg_params);
    s_root_x = secant_solver(egg_wrapper_2x,2,2.75); % WILL CRASH IF DIFFERENCE BETWEEN GUESSES IS >=1
    s_root_x_2 = secant_solver(egg_wrapper_2x,7.75,8.25);
    x_range = sort([s_root_x, s_root_x_2]);
    s_root_y = secant_solver(egg_wrapper_2y,2.25,2.75);
    s_root_y_2 = secant_solver(egg_wrapper_2y,4.25,4.75);
    y_range = sort([s_root_y, s_root_y_2]);
end
%% WRAPPER X
%wrapper function that calls egg_func
%and only returns the x coordinate of the
%point on the perimeter of the egg
%(single output)
function x_out = egg_wrapper_x(s,x0,y0,theta,egg_params)
    [V, G] = egg_func(s,x0,y0,theta,egg_params);
    x_out = G(1);
end
%% WRAPPER Y
%wrapper function that calls egg_func
%and only returns the x coordinate of the
%point on the perimeter of the egg
%(single output)
function y_out = egg_wrapper_y(s,x0,y0,theta,egg_params)
    [V, G] = egg_func(s,x0,y0,theta,egg_params);
    y_out = G(2);
end
%% MAX MIN X
function [x_min, x_max] = egg_min_max(s,x0,y0,theta,egg_params)
    s = 0;
    x_coords = [];
    for n = 1:100
        [V, G] = egg_func(s,x0,y0,theta,egg_params);
        s = s + 0.01;
        x_val = G(1);
        x_coords = [x_coords, x_val];
        x_min = min(x_coords);
        x_max = max(x_coords);
    end
end
%% EGG FUNCTION
%This function generates the parametric curve describing an oval
%INPUTS:
%s: the curve parametr. s is a number from 0 to 1. The curve function has a
% period of 1, so s=.3 and s=1.3 will generate the same output
% s can also be a list (row vector) of numbers
%theta: rotation of the oval. theta is a number from 0 to 2*pi.
% Increasing theta rotates the oval counterclockwise
%x0: horizontal offset of the oval
%y0: vertical offset of the oval
%egg_params: a struct describing the hyperparameters of the oval
% egg_params has three variables, a,b, and c
% without any rotation/translation, the oval satisfies the equation:
% x^2/a^2 + (y^2/b^2)*e^(c*x) = 1
% tweaking a,b,c changes the shape of the oval
%OUTPUTS:
%V: the position of the point on the oval given the inputs
% If s is a single number, then V will have the form of a column vector
% [x_out;y_out] where (x_out,y_out) are the coordinates of the point on
% the oval. If the input t is a list of numbers (a row vector) i.e.:
% s = [s_1,...,s_N]
% then V will be an 2xN matrix:
% [x_1,...,x_N; y_1,...,y_N]
% where (x_i,y_i) correspond to input s_i
%G: the gradient of V taken with respect to s
% If s is a single number, then G will be the column vector [dx/ds; dy/ds]
% If s is the list [s_1,...,s_N], then G will be the 2xN matrix:
% [dx_1/ds_1,...,dx_N/ds_N; dy_1/ds_1,...,dy_N/ds_N]
function [V, G] = egg_func(s,x0,y0,theta,egg_params)
    %unpack the struct
    a=egg_params.a;
    b=egg_params.b;
    c=egg_params.c;
    %compute x (without rotation or translation)
    x = a*cos(2*pi*s);
    %useful intermediate variable
    f = exp(-c*x/2);
    %compute y (without rotation or translation)
    y = b*sin(2*pi*s).*f;
    %compute the derivatives of x and y (without rotation or translation)
    dx = -2*pi*a*sin(2*pi*s);
    df = (-c/2)*f.*dx;
    dy = 2*pi*b*cos(2*pi*s).*f + b*sin(2*pi*s).*df;
    %rotation matrix corresponding to theta
    R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
    %compute position and gradient for rotated + translated oval
    V = R*[x;y]+[x0*ones(1,length(theta));y0*ones(1,length(theta))];
    G = R*[dx;dy];
end