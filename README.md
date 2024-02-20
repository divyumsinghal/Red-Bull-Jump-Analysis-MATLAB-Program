# Red-Bull-Jump-Analysis-MATLAB-Program

This MATLAB program analyzes the Red Bull jump data and performs various computations related to the free fall. The program reads data from a CSV file (`RedBullJumpData.csv`), processes the information, and generates plots to visualize the results.

## Computation Part

### Part 1: Measured Velocity vs Time
The program plots the measured velocity vs time using a red line.

```matlab
figure(1);
h_part1 = plot(t_redbull, v_redbull, 'r-x', 'linewidth', 2.0);
xlabel('Time (secs)');
ylabel('Velocity (m/s)');
title('Measured Velocity vs Time');
grid on;
axis([0 180 0 400]);
set(gca, 'FontSize', 24);
hold on;
```

### Part 2: Linear Freefall Velocity vs Time
The program plots the linear freefall velocity vs time using a black dashed line.

```matlab
h_part2 = plot(t_redbull, v_freefall, 'k--', 'linewidth', 2.0);
```

### Part 3: Atmosphere Entry Time
Calculates the time when the jumper hits the Earth's atmosphere.

```matlab
hit_instant = ... % Calculated time when the jumper hits the atmosphere
fprintf('Mr. B hits the earth''s atmosphere at %f secs after he jumps\n', hit_instant);
```

### Part 4: Numerical Solution Using Euler's Method
Numerically solves the differential equation using Euler's method and plots the result.

```matlab
% Euler's method for numerical solution
% Plot the result using green dashed line with markers
h_part4 = plot(t_euler, v_numerical_1, 'g--o', 'LineWidth', 2.0, 'MarkerSize', 2.5);
```

### Part 5: Percentage Error Calculation
Calculates and prints the percentage error at specific time instants.

```matlab
% Calculate and print percentage error
fprintf('The percentage error at 69 and 180 secs is %1.1f and %3.1f respectively\n', per_error(1), per_error(2));
```

### Part 6: Heun's Method with Variable Drag Coefficient
Applies Heun's method to solve the differential equation with a variable drag coefficient and plots the result.

```matlab
% Heun's method for numerical solution with variable drag coefficient
% Plot the result using black dashed line with markers
h_part6 = plot(t_heun, v_numerical_2, 'k--+', 'LineWidth', 2.0, 'MarkerSize', 6);
```

## Physics Formulas and Methods Explained
- **Linear Freefall Velocity**: \( v_{\text{freefall}} = g \cdot t_{\text{redbull}} \)
- **Atmosphere Entry Calculation**: Detects the time when the measured velocity differs significantly from the linear freefall velocity.
- **Euler's Method**: Numerical method to solve first-order ordinary differential equations.
- **Heun's Method**: Improved numerical method for solving differential equations with a variable drag coefficient.

## Author

Divyum Singhal

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
