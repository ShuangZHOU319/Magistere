clc;
clear all;

%************************************************************* PARAMETERS *********************************************************
d0 = 1;  % reference distance (meter)
PLE = 2; % path loss exponent
target_x = 11.0; % true location x-coordinate (meter)
target_y = 10.0; % true location y-coordinate (meter)
s_x = [5 5 5 8 8 12 12 15 15 15]; % x coordinate of anchor sensors (meter)
s_y = [5 10 15 0 20 0 20 5 10 15]; % y coordinate of anchor sensors (meter)
sigma2 = 1; % noise variance (dB)

% Number of anchors
num_anchors = length(s_x);

%********************************************************** Chapman-Robbins Bound Calculation ********************************************************
% Define perturbations (u vectors)
perturbations = [0.1, 0; -0.1, 0; 0, 0.1; 0, -0.1]; % Example perturbations around the target
num_perturbations = size(perturbations, 1);

% Compute mean RSS values and their derivatives with respect to x and y
for k = 1:num_anchors
    dR_dx(k) = (-10 * PLE / log(10)) * (target_x - s_x(k)) / (norm([target_x; target_y] - [s_x(k); s_y(k)], 2)^2);
    dR_dy(k) = (-10 * PLE / log(10)) * (target_y - s_y(k)) / (norm([target_x; target_y] - [s_x(k); s_y(k)], 2)^2);
end
dR_dx = dR_dx';
dR_dy = dR_dy';

% Compute Chapman-Robbins Bound
chrb_bound = inf;
for i = 1:num_perturbations
    u = perturbations(i, :)';
    numerator = (u' * u)^2;
    denominator = (dR_dx' * u)^2 + (dR_dy' * u)^2;
    chrb_bound_i = numerator / denominator;
    if chrb_bound_i < chrb_bound
        chrb_bound = chrb_bound_i;
    end
end
chrb_bound = sqrt(chrb_bound);

% Display results
disp(['Chapman-Robbins Bound: ', num2str(chrb_bound)]);

%******************************************************************* Plot *****************************************************************
% Create a plot of the anchors and the target
figure;
hold on;
plot(s_x, s_y, 'bo', 'MarkerFaceColor', 'b'); % Plot anchor nodes
plot(target_x, target_y, 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Plot target node
xlabel('X Coordinate (meters)');
ylabel('Y Coordinate (meters)');
title('Anchor Nodes and Target Node');
text(target_x, target_y, ' Target', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'red', 'FontSize', 10);
for k = 1:num_anchors
    text(s_x(k), s_y(k), sprintf(' Anchor %d', k), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'blue', 'FontSize', 8);
end
grid on;
hold off;
