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

%********************************************************** Barankin Bound Calculation ********************************************************
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

% Fisher Information Matrix (FIM)
FIM = zeros(2, 2);
FIM(1, 1) = (dR_dx' * dR_dx) / sigma2;
FIM(2, 2) = (dR_dy' * dR_dy) / sigma2;
FIM(1, 2) = (dR_dx' * dR_dy) / sigma2;
FIM(2, 1) = FIM(1, 2);

% Compute Barankin Bound
barankin_bound = inf;
for i = 1:num_perturbations
    u = perturbations(i, :)';
    u_tilde = FIM * u;
    barankin_bound_i = (u' * u) / (u_tilde' * u_tilde);
    if barankin_bound_i < barankin_bound
        barankin_bound = barankin_bound_i;
    end
end
barankin_bound = sqrt(barankin_bound);

% Display results
disp(['Barankin Bound: ', num2str(barankin_bound)]);

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
