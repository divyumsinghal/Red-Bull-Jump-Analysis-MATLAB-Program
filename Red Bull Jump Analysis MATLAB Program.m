%Red Bull Data plotted.
close all
clear 

% Reading data from file
jumpdata = csvread('RedBullJumpData.csv');

% First column is time in secs
t_redbull = jumpdata(:,1);

% Second column is velocity

v_redbull = jumpdata(:,2);

% Third column is theoretical terminal velocity at that time
terminal_velocity = jumpdata(:,3); % You need to use this in the last question

%Storing number of rows in a variable
N_timestamps = length(t_redbull);


%Calculated freefall velocity vector here
g = 9.81;
v_freefall = g * t_redbull;

% Part 1

%Drawing the red line for the measured velocity and the time
figure(1);
h_part1 = plot(t_redbull, v_redbull, 'r-x', 'linewidth', 2.0);
xlabel('Time (secs)');
ylabel('Velocity (m/s)');
title('Measured Velocity vs Time');

hold on;

% Part 2
%Drawing the black dashed line for the linear freefall velocity and the time
h_part2 = plot(t_redbull, v_freefall, 'k--', 'linewidth', 2.0); 
shg;

%Some nessacary graph instuctions
grid on;
axis([0 180 0 400]);
set(gca, 'FontSize', 24);
hold on;

 
% Part  3
% Calculating when he hits the atmosphere
%finding the index in the arrays where the difference PERCENTAGE (hence the need for 0.02 * v_redbull(2:end), instead of just 0.02) between the values is too large
%this will not give the exact instant as the values in csv file are not continous, but rather in increments. This gives the first increment to satisfy the conditions.

N = length(t_redbull);

for i = 2:1:N
    if (abs(v_redbull(i) - v_freefall(i)) / v_redbull(i) >= 0.02)
        hit_instant = t_redbull(i);
        break
    end
end

fprintf('Mr. B hits the earth''s atmoshpere at %f secs after he jumps\n', hit_instant);

% Part 4
% Now starting from the velocity and time = 56 secs ; index 16 stored in variable start

start = find(t_redbull == 56);
drag_constant = 3/60;

% Starting from this time instant, calculating the velocity required and now the time vector for plotting the euler's solution
% t_euler is an array from 56 sec to the end of t_redbull

t_euler = t_redbull(start:end);
N = length(t_euler);

% Initializing an array to store numerical results with N elements
v_numerical_1 = zeros(N, 1);  

% Setting the initial value of the array using a function v_redbull at the start position
v_numerical_1(1) = v_redbull(start);

for i = 1:1:N-1

    % Calculating the derivative dv/dt using a simple model: g - drag_constant * v_numerical_1(i-1)
    dvdt = g - drag_constant * v_numerical_1(i);
    
    % Using Euler's method to update the velocity for the next time step
    v_numerical_1(i+1) = v_numerical_1(i) + dvdt * (t_euler(i+1) - t_euler(i));

end

% Ploting the first numerical solution v/s time using the dashed green line with (+) markers
h_part4 = plot(t_euler, v_numerical_1, 'g--o', 'LineWidth', 2.0, 'MarkerSize', 2.5);
shg;

% Part 5 

%get the indices of each point in t_euler & t_redbull
% technically, once you have t_euler, you can get t_redbull due to their direct relation, but still

index_t_euler_first = find(t_euler == 69); 
index_t_euler_last = find(t_euler == 180); 
index_t_redbull_first = find(t_redbull == 69); 
index_t_redbull_last = find(t_redbull == 180); 

% Calculate the percentage error as required & then print them

per_error = [100 * abs(v_redbull(index_t_redbull_first) - v_numerical_1(index_t_euler_first)) / v_redbull(index_t_redbull_first), 100 * abs(v_redbull(index_t_redbull_last) - v_numerical_1(index_t_euler_last)) / v_redbull(index_t_redbull_last)];

fprintf('The percentage error at 69 and 180 secs is %1.1f and\n', per_error(1));
fprintf('%3.1f respectively\n', per_error(2));


% Part 6 
% You'll need to repeat your euler loop here again but this time
% update the drag constant at every timestamp and change the update
% calculation to allow for the new v^2(t) term
% A hint here that now you have to calculate the velocity using the new
% differental equation
% constant .. put it in v_numerical_2
% This is the handle plot for part 6. You have to plot the right stuff not
% this stuff.
% Note that the plot linewidth and colour are wrong. Fix it.
% Solving the new differential equation using Heun's solution

heun_end = find(t_redbull == 100);
t_heun = t_redbull(start:heun_end);

N = length(t_heun);

v_numerical_2 = zeros(N, 1);
v_numerical_2(1) = v_redbull(start);

for i = 1:1:N-1

    % Calculate c(t)/m using the measured terminal velocity
    % by stands for divided by
    % terminal_velocity are shofted by 1 in the original dataset

    ct_by_m = g / terminal_velocity(start + i )^2;

    % Heun's method to update the velocity for the next time step
    % Technically I am using 2nd order Range Kutta with a1 = 1/2; a2 = 1/2 ; p1 = q1 = 1 as it is the same as Heun's method and very easy to implement in a for loop

    %k1 = f(xi, yi)

    k1 = g - ct_by_m * v_numerical_2(i)^2;

    % Updating my ct_by_m  for the k2 values

    ct_by_m = g / terminal_velocity(start + i + 1)^2;

    %k2 = f(xi + h, yi + k1h)

    k2 = g - ct_by_m * (v_numerical_2(i) + k1 * (t_heun(i+1) - t_heun(i)))^2;

    dvdt = (1/2) * k1 + (1/2) * k2;

    % yi+1 = yi + (a1k1 + a2k2)h

    v_numerical_2(i+1) = v_numerical_2(i) + dvdt * (t_heun(i+1) - t_heun(i));

end


% Plot using dashed black line with (+) markers
h_part6 = plot(t_heun, v_numerical_2, 'k--+', 'LineWidth', 2.0, 'MarkerSize', 6);
shg;

% Calculate the percentage error at t = 100 secs
index_t_heun_100s = find(t_heun == 100);
index_t_redbull_100s = find(t_redbull == 100);

est_error = 100 * abs(v_redbull(index_t_redbull_100s) - v_numerical_2(index_t_heun_100s)) / v_redbull(index_t_redbull_100s);

fprintf('The error at t = 100 secs using my estimated drag information is %f\n', est_error);


% Printing all values for debegging and testing

% v_freefall
% hit_instant
% per_error
% v_numerical_1
% est_error
% v_numerical_2
% h_part1
% h_part2
% h_part4
% h_part6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE THIS IS TO MAKE SURE YOU HAVE USED THE
% VARIABLES THAT WE ASKED FOR
% Check for existence of variables
if (~exist('v_freefall', 'var'))
  error('The variable v_freefall does not exist.')
end;
if (~exist('hit_instant', 'var'))
  error('The variable hit_instant does not exist.')
end;
if (~exist('per_error', 'var'))
  error('The variable per_error does not exist.')
end;
if (exist('per_error', 'var'))
  l = size(per_error);
  if ( sum(l - [1 2]) ~= 0)
    error('per_error is not a 2 element vector. Please make it so.')
  end;
end;
if (~exist('v_numerical_1', 'var'))
  error('The variable v_numerical_1 does not exist.')
end;  
if (~exist('est_error', 'var'))
  error('The variable est_error does not exist.')
end;  
if (~exist('h_part1', 'var'))
  error('The plot handle h_part11 is missing. Please create it as instructed.')
end;
if (~exist('h_part2', 'var'))
  error('The plot handle h_part11 is missing. Please create it as instructed.')
end;
if (~exist('h_part4', 'var'))
  error('The plot handle h_part11 is missing. Please create it as instructed.')
end;
if (~exist('h_part6', 'var'))
  error('The plot handle h_part11 is missing. Please create it as instructed.')
end;


