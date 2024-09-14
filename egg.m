egg_params = struct();
egg_params.a = 3; egg_params.b = 2; egg_params.c = .15;
%specify the position and orientation of the egg
x0 = 10; y0 = 10; theta = pi/6;
figure;
%compute the perimeter of the egg
% hold on; axis equal; axis square
% axis([0,10,0,10])
% [V_list, G_list] = egg_func(linspace(0,1,100),x0,y0,theta,egg_params);
% %plot the perimeter of the egg
% plot(V_list(1,:),V_list(2,:),'k');
% hold on
% [x_min, x_max, y_min, y_max] = bounding_box(x0, y0, theta, egg_params)
% plot([xmin, xmax, xmax, xmin, xmin], [ymin ymin, ymax, ymax, ymin])
% hold off
animation_example(x0,y0,theta,egg_params,@egg_trajectory01,30,0)
%%
traj_fun = @egg_trajectory01
y_ground = 2
x_wall = 30
[t_ground,t_wall] = collision_func(@egg_trajectory01, egg_params, y_ground, x_wall)

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
%% LEFT RIGHT BOUNDS
function [xmin, xmax] = LR_bounding(x0, y0, theta, egg_params)
    egg_wrapper_2x = @(s) egg_wrapper_x(s,x0,y0,theta,egg_params);
    s_guess_list = [0, 1/3, 2/3];
    x_list = zeros(length(s_guess_list), 1);
    for n=1:length(s_guess_list)
        s_guess = s_guess_list(n);
        s_root_x = secant_solver(egg_wrapper_2x, s_guess, s_guess+1e-5);
        [V,G] = egg_func(s_root_x, x0, y0, theta, egg_params);
        x_list(n) = V(1);
    end
    xmin = min(x_list);
    xmax = max(x_list);
end
%% TOP BOTTOM BOUNDS
function [ymin, ymax] = TB_bounding(x0, y0, theta, egg_params)
    egg_wrapper_2y = @(s) egg_wrapper_y(s,x0,y0,theta,egg_params);
    s_guess_list = [0, 1/3, 2/3];
    y_list = zeros(length(s_guess_list), 1);
    for n=1:length(s_guess_list)
        s_guess = s_guess_list(n);
        s_root_y = secant_solver(egg_wrapper_2y, s_guess, s_guess+1e-5);
        [V,G] = egg_func(s_root_y, x0, y0, theta, egg_params);
        y_list(n) = V(2);
    end
    ymin = min(y_list);
    ymax = max(y_list);
end

%% BOUNDING BOX
function [x_min, x_max, y_min, y_max] = bounding_box(x0, y0, theta, egg_params)
    [x_min, x_max] = LR_bounding(x0, y0, theta, egg_params);
    [y_min, y_max] = TB_bounding(x0, y0, theta, egg_params);
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
%% PARABOLIC TRAJECTORY
function [x0,y0,theta] = egg_trajectory01(t)
    x0 = 7*t + 8;
    y0 = -6*t.^2 + 20*t + 6;
    theta = 0.05*t;
end
%% BOUNDING WRAPPER X
function x_traj = x_bounding_traj(t, traj_fun, egg_params)
    [x0,y0,theta] = traj_fun(t);
    [x_range, ~] = bounding_box(x0, y0, theta, egg_params);
    x_traj = max(x_range);
end
%% BOUNDING WRAPPER Y
function y_traj = y_bounding_traj(t, traj_fun, egg_params)
    [x0,y0,theta] = traj_fun(t);
    [~ ,y_range] = bounding_box(x0, y0, theta, egg_params);
    y_traj = y_range(1);
end
%% COLLISION
%Function that computes the collision time for a thrown egg
%INPUTS:
%traj_fun: a function that describes the [x,y,theta] trajectory
% of the egg (takes time t as input)
%egg_params: a struct describing the hyperparameters of the oval
%y_ground: height of the ground
%x_wall: position of the wall
%OUTPUTS:
%t_ground: time that the egg would hit the ground
%t_wall: time that the egg would hit the wall

function [t_ground,t_wall] = collision_func(traj_fun, egg_params, y_ground, x_wall)
    x_traj_wrapper_2 = @(t) x_bounding_traj(t, traj_fun, egg_params) - x_wall;
    y_traj_wrapper_2 = @(t) y_bounding_traj(t, traj_fun, egg_params) - y_ground;
%     t_ground = bisection_solver(y_traj_wrapper_2, 0, 5)
%     t_wall = 2
    if bisection_solver(x_traj_wrapper_2, 0, 10) ~= 0
        t_ground = bisection_solver(y_traj_wrapper_2, 0, 5);
        t_wall = "NULL";
    elseif bisection_solver(y_traj_wrapper, 0, 5) ~= 0
        t_ground = "NULL";
        t_wall = bisection_solver(x_traj_wrapper_2, 0, 10);
    else
        t_ground = bisection_solver(y_traj_wrapper_2, 0, 5);
        t_wall = bisection_solver(x_traj_wrapper_2, 0, 10);
    end
end
%% ANIMATION
%Short example demonstrating how to create a MATLAB animation
%In this case, a square moving along an elliptical path
function animation_example(x0,y0,theta,egg_params,traj_fun,x_wall,y_ground)

    hold on; axis equal; axis square
    axis([0,40,0,40])
    [V_list, G_list] = egg_func(linspace(0,1,100),x0,y0,theta,egg_params);
    x_coords = V_list(1,:);
    y_coords = V_list(2,:);
    egg_plot = plot(x_coords, y_coords,'k');
    hold on
    [x_min, x_max, y_min, y_max] = bounding_box(x0, y0, theta, egg_params);
    [t_ground,t_wall] = collision_func(traj_fun, egg_params, y_ground, x_wall)
    for t=0:.001:10
        [V_list, G_list] = egg_func(linspace(0,1,100),x0,y0,theta,egg_params);
        [x_shift, y_shift, theta_shift] = traj_fun(t);
        x_coords = V_list(1,:);
        y_coords = V_list(2,:);
        position_x = x_shift;
        position_y = y_shift;
        theta = theta + theta_shift;
        %compute positions of square vertices (in world frame)
        x_plot = (x_coords+position_x);
        y_plot = (y_coords+position_y);
        yline(y_ground)
        xline(x_wall)
        %update the coordinates of the square plot
        set(egg_plot,'xdata',x_plot,'ydata',y_plot);
        %update the actual plotting window
        drawnow;
        if t == t_ground || t == t_wall
            pause(2)
            break
        end
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