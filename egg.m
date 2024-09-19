%% FULL EGG ANIMATION
function egg()

    % specify parameters to make egg
    egg_params = struct();
    egg_params.a = 3; egg_params.b = 2; egg_params.c = .15;

    % animate collision of egg and save video
    animation(egg_params,@egg_trajectory01,10,0)
end
%% LEFT RIGHT BOUNDS
function [xmin, xmax, s_root_x] = LR_bounding(x0, y0, theta, egg_params)

    % get one point on the egg curve
    egg_wrapper_2x = @(s) egg_wrapper_x(s,x0,y0,theta,egg_params);
    s_guess_list = [0, 1/3, 2/3];
    x_list = zeros(length(s_guess_list), 1);

    % find x-coordinates along the edge of the egg
    for n=1:length(s_guess_list)
        s_guess = s_guess_list(n);
        s_root_x = secant_solver(egg_wrapper_2x, s_guess, s_guess+1e-5);
        [V,G] = egg_func(s_root_x, x0, y0, theta, egg_params);
        x_list(n) = V(1);
    end

    % determine the upper and lower x bounds
    xmin = min(x_list);
    xmax = max(x_list);
end
%% TOP BOTTOM BOUNDS
function [ymin, ymax, s_root_y] = TB_bounding(x0, y0, theta, egg_params)

    % get one point on the egg curve
    egg_wrapper_2y = @(s) egg_wrapper_y(s,x0,y0,theta,egg_params);
    s_guess_list = [0, 1/3, 2/3];
    y_list = zeros(length(s_guess_list), 1);

    % find y-coordinates along the edge of the egg
    for n=1:length(s_guess_list)
        s_guess = s_guess_list(n);
        s_root_y = secant_solver(egg_wrapper_2y, s_guess, s_guess+1e-5);
        [V,G] = egg_func(s_root_y, x0, y0, theta, egg_params);
        y_list(n) = V(2);
    end

    % determine the upper and lower y bounds
    ymin = min(y_list);
    ymax = max(y_list);
end

%% BOUNDING BOX
function [x_min, x_max, y_min, y_max] = bounding_box(x0, y0, theta, egg_params)

    % find min and max x and y values
    [x_min, x_max] = LR_bounding(x0, y0, theta, egg_params);
    [y_min, y_max] = TB_bounding(x0, y0, theta, egg_params);
end

%% WRAPPER X
function x_out = egg_wrapper_x(s,x0,y0,theta,egg_params)

    [V, G] = egg_func(s,x0,y0,theta,egg_params);

    % x-coordinate of single point on egg curve
    x_out = G(1);
end
%% WRAPPER Y
function y_out = egg_wrapper_y(s,x0,y0,theta,egg_params)

    [V, G] = egg_func(s,x0,y0,theta,egg_params);

    % y-coordinate of single point on egg curve
    y_out = G(2);
end
%% PARABOLIC TRAJECTORY
function [x0,y0,theta] = egg_trajectory01(t)

    % functions to determine the trajectory of the egg
    x0 = 10*t -10;
    y0 = -6*t.^2 + 20*t + 3;
    theta = 2*t;
end
%% BOUNDING WRAPPER X
function x_traj = x_bounding_traj(t, traj_fun, egg_params)
    [x0,y0,theta] = traj_fun(t);
    [x_min, x_max, y_min, y_max] = bounding_box(x0, y0, theta, egg_params);

    % x-coordinate of furthest (rightmost) point on the bounding box at given time
    x_traj = x_max;
end
%% BOUNDING WRAPPER Y
function y_traj = y_bounding_traj(t, traj_fun, egg_params)
    [x0,y0,theta] = traj_fun(t);
    [x_min, x_max, y_min, y_max] = bounding_box(x0, y0, theta, egg_params);

    % y-coordinate of lowest point on the bounding box at given time
    y_traj = y_min;
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

function [t_ground,t_wall] = collision_func(traj_fun, egg_params, x_wall, y_ground)

    % find where distance between the wall/ground and bounding box is zero
    x_traj_wrapper_2 = @(t) x_bounding_traj(t, traj_fun, egg_params) - x_wall;
    y_traj_wrapper_2 = @(t) y_bounding_traj(t, traj_fun, egg_params) - y_ground;

    % find the time to collide with the ground/wall using bisection method
    t_ground = bisection_solver(y_traj_wrapper_2, 0, 10);
    t_wall = bisection_solver(x_traj_wrapper_2, 0, 10); 
end
%% ANIMATION
function animation(egg_params,traj_fun,x_wall,y_ground)
    
    % Download path for egg animation
    mypath1 = 'C:\Users\ldao\Downloads\';
    fname='egg_animation.avi';
    input_fname = [mypath1,fname];

    % create a videowriter, which will write frames to the animation file
    writerObj = VideoWriter(input_fname);

    % must call open before writing any frames
    open(writerObj);
    
    % make egg and plot
    [x0,y0,theta] = traj_fun(0);
    fig1 = figure(1);
    hold on; axis equal; axis square
    axis([-20,40,-20,40])
    [V_list, G_list] = egg_func(linspace(0,1,100),x0,y0,theta,egg_params);
    x_coords = V_list(1,:);
    y_coords = V_list(2,:);
    egg_plot = plot(x_coords, y_coords,'k');
    hold on

    % determine the time the egg will collide
    [t_ground,t_wall] = collision_func(traj_fun, egg_params, x_wall, y_ground);
    stop = min([t_wall,t_ground]);

    % plot ground and wall
    yline(y_ground);
    xline(x_wall);

    % plot egg moving along trajectory until collision
    for t=0:.01:3

        % trajectory path
        [x_shift, y_shift, theta_shift] = traj_fun(min(t, stop));
        
        [V_list, G_list] = egg_func(linspace(0,1,100),x_shift,y_shift,theta_shift,egg_params);

        %[x_min, x_max, y_min, y_max] = bounding_box(x0, y0, theta, egg_params);
        
        % find coordinates of egg
        x_coords = V_list(1,:);
        y_coords = V_list(2,:);
    
        % plot egg
        set(egg_plot,'xdata',x_coords,'ydata',y_coords);
        drawnow;
        hold on

        % plot point where egg collides
        if t >= stop
            % collision with ground
            if t_ground < t_wall
                y_stop = y_bounding_traj(stop, traj_fun, egg_params);
                [x_stop, y_pt, theta] = traj_fun(stop);
                plot(x_stop, y_stop, 'o', 'MarkerFaceColor', 'r')
            % collision with wall
            elseif t_wall < t_ground
                x_stop = x_bounding_traj(stop, traj_fun, egg_params);
                [x_pt, y_stop, theta] = traj_fun(stop);
                plot(x_stop, y_stop, 'o', 'MarkerFaceColor', 'r')
            end
        end

        current_frame = getframe(fig1);

        % write the frame to the video
        writeVideo(writerObj,current_frame);

    end
    
    % must call close after all frames are written
    close(writerObj);

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