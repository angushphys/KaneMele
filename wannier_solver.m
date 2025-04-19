%%
% Plot Wannier result
% This program requires two inputs:
% - KPOINTS: path for band structure calculation
% - wan_basis.mat: Wannier hopping parameters

clear;

%%

% Set the Fermi level
E_f = 0;

% Set the energy window (relative to Fermi level)
E_min = -5;
E_max =  5;

% Number of k-points for each segment between high-symmetry points
n_k_each_line = 201;

% Names of high-symmetry points (used in band structure plot)
x_tic_name = {'G' 'X' 'K' 'G' 'L' 'W' 'X' 'K' 'W'};

% Label formatting options for x-axis
% 0 = none
% 1 = bar
% 2 = tilde
% 3 = hat
% 4 = check
index_add_bar = 0;

%%

% Load Wannier basis and hopping parameters
load wan_basis
lattice_a = wan_basis.lattice_a;
lattice_b = wan_basis.lattice_b;
hopping = wan_basis.hopping;

% % (Commented out) Transformation from direct to Cartesian coordinates
% % k_cart_all = zeros(3, n_k);
% % for i_k = 1 : n_k
% %     k_cart_all(:, i_k) = lattice_b' * k_direct_all(:, i_k);
% % end

% % (Commented out) Transformation from Cartesian to direct coordinates
% % k_diret_all = zeros(3, n_k);
% % for i_k = 1 : n_k
% %     k_direct_all(:, i_k) = lattice_b \ k_cart_all(:, i_k);
% % end

%%

% Read KPOINTS file in line-mode to extract k-path for band structure
clear high_sym_points;
f_kp = fopen('KPOINTS');
if f_kp <= 0
    error('==  Could not open KPOINTS file!  ==');
end

% Skip header
for i_line = 1:3
    fgets(f_kp);
end

% Read coordinate mode (Reciprocal/Cartesian)
kpoint_coordinate_mode_char = fscanf(f_kp, '%c', 1);

if strcmp(kpoint_coordinate_mode_char, 'R') || ...
   strcmp(kpoint_coordinate_mode_char, 'r') || ...
   strcmp(kpoint_coordinate_mode_char, 'D') || ...
   strcmp(kpoint_coordinate_mode_char, 'd')
    
    k_point_coordinate_mode = 1;  % Direct or Reciprocal coordinates

elseif strcmp(kpoint_coordinate_mode_char, 'C') || ...
       strcmp(kpoint_coordinate_mode_char, 'c')

    k_point_coordinate_mode = 0;  % Cartesian coordinates
    
else
    error('== ERROR: Unknown coordinate mode in KPOINTS ==');
end

% Read high-symmetry points
fgets(f_kp);
high_sym_points(1:3, 1) = fscanf(f_kp, '%f', 3)';
fgets(f_kp);
n_high_sym_points = 1;
while 1
    [sym_point, scan_count] = fscanf(f_kp, '%f', 3);
    if scan_count == 0
        break;
    end

    high_sym_points(:, n_high_sym_points+1) = sym_point;

    n_high_sym_points = n_high_sym_points + 1;
    fgets(f_kp);

    fscanf(f_kp, '%f', 1);
    fgets(f_kp);
end
fclose(f_kp);
clear f_kp kp_mod scan_count;

n_high_sym_points = size(high_sym_points, 2);

% Convert Cartesian to direct coordinates if needed
if k_point_coordinate_mode == 0
    for i_k = 1 : n_high_sym_points
        high_sym_points(:, i_k) = lattice_b \ high_sym_points(i_k, :);
    end
end

% Total number of interpolated k-points
n_k = n_k_each_line * (n_high_sym_points - 1);
k_point_all = zeros(3, n_k);

% Fill in default names for any extra high-symmetry points
if size(x_tic_name, 2) < n_high_sym_points
    for i_tic_name = size(x_tic_name, 2)+1 : n_high_sym_points
        x_tic_name{i_tic_name} = 'X';
    end
end

% Linearly interpolate k-points between high-symmetry points
for i_k = 1:(n_high_sym_points-1)
    for i_direction = 1 : 3
        k_begin = n_k_each_line*(i_k - 1) + 1;
        k_end = n_k_each_line*i_k;
        
        k_point_all(i_direction, k_begin : k_end) = linspace( ...
            high_sym_points(i_direction, i_k), ...
            high_sym_points(i_direction, i_k+1), ...
            n_k_each_line);
    end
end

% Compute k-path distances
k_axis = zeros(n_k, 1);
k_axis(1) = 0;
n_x_tic = 1;
x_tic = zeros(2,1);

for i_k = 2 : n_k
    distance = norm(...
        lattice_b' * (k_point_all(:, i_k) - k_point_all(:, i_k-1)) );

    if distance == 0.0
        distance = eps*100;

        n_x_tic = n_x_tic + 1;
        x_tic(n_x_tic) = k_axis(i_k-1);
    end
    k_axis(i_k) = k_axis(i_k-1) + distance;
end
n_x_tic = n_x_tic + 1;
x_tic(n_x_tic) = k_axis(i_k);

clear distance;

%%

% Main part: calculate Hamiltonians and eigenvalues
H_k = zeros(wan_basis.n_band, wan_basis.n_band);
H_eig_0 = zeros(n_k, wan_basis.n_band);

fprintf('  Begin eigenvalue calculations...\n');

% Loop over each k-point to compute eigenvalues
for i_k = 1 : n_k
    k_direct = k_point_all(:, i_k);
    
    % Construct Hamiltonian using hopping parameters:
    % Hk(a, b) = sum_R{ t(a, b) * exp(i * 2π * k · R) }
    H_k = full(sparse( ...
            hopping.orbit_0, hopping.orbit_R,...
            exp(2i * pi * hopping.R * k_direct) .* hopping.t, ...
            wan_basis.n_band, wan_basis.n_band));

    H_k = (H_k + H_k') / 2;  % Make Hamiltonian Hermitian

    H_eig_0(i_k, :) = eig(H_k);  % Diagonalize H to get eigenvalues

    if mod(i_k, n_k_each_line) == 0
        fprintf('  Finished calculating %2d segment(s)...\n', i_k/n_k_each_line);
    end

    if n_k_each_line > 300 && mod(i_k, 200) == 0
        fprintf('  ---> Finished %4d k-points...\n', i_k);
    end
end

fprintf('  Eigenvalue calculation completed.\n\n');

%%

% Format tick labels for plotting
% G -> Gamma, Sg -> Sigma, add bar/tilde/hat if specified

for i_x_tic = 1 : n_x_tic
    tic_length = length(x_tic_name{i_x_tic});
    
    if strcmp(x_tic_name{i_x_tic}(1), 'G')
        x_tic_name{i_x_tic} = sprintf('\\Gamma%s',x_tic_name{i_x_tic}(2:tic_length));
    elseif tic_length >= 2 && strcmp(x_tic_name{i_x_tic}(1:2), 'Sg')
        x_tic_name{i_x_tic} = sprintf('\\Sigma%s',x_tic_name{i_x_tic}(3:tic_length));
    end
    tic_length = length(x_tic_name{i_x_tic});
    
    index_sub = 0;
    x_tic_sub = '';
    if tic_length > 2 && x_tic_name{i_x_tic}(tic_length-1) == '_'
        index_sub = 1;
        x_tic_sub = x_tic_name{i_x_tic}(tic_length);
        x_tic_name{i_x_tic} = x_tic_name{i_x_tic}(1:tic_length-2);
    elseif tic_length > 1 && strcmp(x_tic_name{i_x_tic}(tic_length), '''')
        index_sub = 2;
        x_tic_name{i_x_tic} = x_tic_name{i_x_tic}(1:tic_length-1);
        fprintf('2\n');
    end
    
    if index_add_bar == 0
        x_tic_name{i_x_tic} = sprintf('%s', x_tic_name{i_x_tic});
    elseif index_add_bar == 1
        x_tic_name{i_x_tic} = sprintf('\\bar{%s}', x_tic_name{i_x_tic});
    elseif index_add_bar == 2
        x_tic_name{i_x_tic} = sprintf('\\tilde{%s}', x_tic_name{i_x_tic});
    elseif index_add_bar == 3
        x_tic_name{i_x_tic} = sprintf('\\hat{%s}', x_tic_name{i_x_tic});
    elseif index_add_bar == 4
        x_tic_name{i_x_tic} = sprintf('\\check{%s}', x_tic_name{i_x_tic});
    end
    
    if index_sub == 1
        x_tic_name{i_x_tic} = sprintf('%s_%s', x_tic_name{i_x_tic}, x_tic_sub);
    elseif index_sub == 2
        x_tic_name{i_x_tic} = sprintf('%s''', x_tic_name{i_x_tic});
    end
    x_tic_name{i_x_tic} = sprintf('$$\\mathrm{%s}$$', x_tic_name{i_x_tic});
end

%%

% Plot band structure
H_eig = H_eig_0 - E_f;

figure('Position', [20, 60, 560, 600]);
subplot('Position', [0.16, 0.08, 0.80, 0.85]);

for i_x_tic = 2 : (n_x_tic-1)
    plot([x_tic(i_x_tic), x_tic(i_x_tic)], ...
        [E_min-2.0, E_max+2.0], ...
        'Color', [0.7, 0.7, 0.7], ...
        'LineWidth', 1);
    hold on;
end
plot([x_tic(1), x_tic(n_x_tic)], ...
    [0.0, 0.0], ...
    '--',...
    'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(k_axis, H_eig, '-b', 'LineWidth', 2);
axis([-inf, Inf, E_min, E_max]);

ylhand = get(gca, 'ylabel');
set(ylhand, 'string', 'Energy (eV)', ...
    'fontsize', 20, ...
    'color', 'k');
set(gca, ...
    'XTick', x_tic, ...
    'XTickLabel', x_tic_name, ...
    'LineWidth', 1, ...
    'FontSize', 18, ...
    'TickLabelInterpreter', 'latex', ...
    'FontName', 'Times New Roman' , ...
    'TickDir', 'out', ...
    'TickLength', 1.5*get(gca,'TickLength') );
