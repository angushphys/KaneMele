% For MATLAB 2016b or later
% Kane-Mele model
% C. L. Kane and E. J. Mele, Phys. Rev. Lett. 95, 226801 (2005)

hopping_on_site_A = 0;
hopping_on_site_B = 0;
hopping_nn_t = -3;

hopping_nnn_kane_mele = -0.05;
hopping_nnn_t = 0;

%%
% Graphene structure parameters
graphene_bond_length = 1.42;

% z-axis length in Angstrom
lattice_a_3 = 20.0;

% Lattice constant in Angstrom
lattice_const = graphene_bond_length * sqrt(3);

% Lattice vectors (in Angstrom units)
lattice_a_0 = [ ...
    1.0000   0.0000
   -0.5000   sqrt(3)/2  ];
lattice_a = zeros(3);
lattice_a(1:2, 1:2) = lattice_const .* lattice_a_0;
lattice_a(3, 3) = lattice_a_3;

% Reciprocal lattice vectors (in 1/Angstrom)
lattice_b = inv(lattice_a)' * 2 * pi;

% Number of ions per unit cell
n_ion_unit_cell = 2;

% Ion positions in fractional (direct) coordinates
ion_position_all_xy = [ ...
    1/3  2/3
    2/3  1/3 ];
ion_position_all = zeros(n_ion_unit_cell, 3);
ion_position_all(:, 1:2) = ion_position_all_xy;

% Ion positions in Cartesian coordinates (Angstrom units)
ion_cartesian_all = zeros(n_ion_unit_cell, 3);
for i_ion = 1 : n_ion_unit_cell
    ion_cartesian_all(i_ion, :) = ...
        ion_position_all(i_ion, :) * lattice_a;
end

%%
% Write POSCAR file (VASP format) for visualization or verification
is_write_poscar = 0;
if is_write_poscar
    f_poscar_name = 'POSCAR';
    f_poscar = fopen(f_poscar_name, 'w');
    fprintf(f_poscar, 'Graphene\n');
    fprintf(f_poscar, '%13.9f\n', lattice_const);
    for i_lattice_direction = 1 : 3
        fprintf(f_poscar, ' %13.9f %13.9f %13.9f\n', ...
            lattice_a(i_lattice_direction, 1) ./ lattice_const, ...
            lattice_a(i_lattice_direction, 2) ./ lattice_const, ...
            lattice_a(i_lattice_direction, 3) ./ lattice_const);
    end
    fprintf(f_poscar, '    C\n');
    fprintf(f_poscar, '   %2d\n', n_ion_unit_cell);
    fprintf(f_poscar, ' Direct\n');
    for i_ion = 1 : n_ion_unit_cell
        fprintf(f_poscar, '%13.9f %13.9f %13.9f\n', ...
            ion_position_all(i_ion, 1), ...
            ion_position_all(i_ion, 2), ...
            ion_position_all(i_ion, 3));
    end
    fclose(f_poscar);
end

%%
% Define the region in real space to search for hopping terms
R_search_x_range = 1;
R_search_y_range = 1;
R_search_x_max =  R_search_x_range;
R_search_y_max =  R_search_y_range;
R_search_x_min = -R_search_x_range;
R_search_y_min = -R_search_y_range;

% Generate all lattice translation vectors R within the search region
n_R_search_x = R_search_x_range * 2 + 1;
n_R_search_y = R_search_y_range * 2 + 1;
n_R_search = n_R_search_x * n_R_search_y;
R_search_all = zeros(n_R_search, 3);
i_R = 1;
for i_R_x = R_search_x_min : R_search_x_max
    for i_R_y = R_search_y_min : R_search_y_max
        R_search_all(i_R, :) = [i_R_x i_R_y 0];
        i_R = i_R + 1;
    end
end

%%
% Generate a list of all atomic positions in the search region
% The 3rd dimension indicates the atom type: A = 1, B = 2
ion_position_search_all = zeros(n_R_search, 3, n_ion_unit_cell);
for i_ion = 1 : n_ion_unit_cell
    ion_position_search_all(:, :, i_ion) = ...
        R_search_all + ion_position_all(i_ion, :);
end

% Convert to Cartesian coordinates
ion_cartesian_search_all = zeros(n_R_search, 3, n_ion_unit_cell);
for i_ion = 1 : n_ion_unit_cell
    for i_R = 1 : n_R_search
        ion_cartesian_search_all(i_R, :, i_ion) = ...
            ion_position_search_all(i_R, :, i_ion) * lattice_a;
    end
end

%%
% Identify nearest-neighbor (NN) and next-nearest-neighbor (NNN) atoms based on distance
nn_distance = graphene_bond_length;
nnn_distance = graphene_bond_length * sqrt(3);

% Search for NN (A ↔ B)
ion_nn_R_each_site = cell(n_ion_unit_cell, 1);
ion_nn_cartesian_each_site = cell(n_ion_unit_cell, 1);
for i_site = 1 : n_ion_unit_cell
    if i_site == 1
        i_site_R = 2;
    else
        i_site_R = 1;
    end
    r_between_ion = ...
        ion_cartesian_search_all(:, :, i_site_R) ...
        - ion_cartesian_all(i_site, :);
    distance_between_ion = sqrt(sum(r_between_ion.^2, 2));
    hopping_nn_index = abs(distance_between_ion - nn_distance) < 1e-7;
    ion_nn_R_each_site{i_site} = R_search_all(hopping_nn_index, :);
    ion_nn_cartesian_each_site{i_site} = ...
        ion_cartesian_search_all(hopping_nn_index, :, i_site_R);
end

% Search for NNN (A ↔ A, B ↔ B)
ion_nnn_R_each_site = cell(n_ion_unit_cell, 1);
ion_nnn_cartesian_each_site = cell(n_ion_unit_cell, 1);
for i_site = 1 : n_ion_unit_cell
    i_site_R = i_site;
    r_between_ion = ...
        ion_cartesian_search_all(:, :, i_site_R) - ion_cartesian_all(i_site, :);
    distance_between_ion = sqrt(sum(r_between_ion.^2, 2));
    hopping_nnn_index = abs(distance_between_ion - nnn_distance) < 1e-7;
    ion_nnn_R_each_site{i_site} = R_search_all(hopping_nnn_index, :);
    ion_nnn_cartesian_each_site{i_site} = ...
        ion_cartesian_search_all(hopping_nnn_index, :, i_site_R);
end

%%
% Debug visualization (optional)
is_debug = 0;
if is_debug == 1
    figure
    for i_ion = 1
        scatter( ...
            ion_nn_cartesian_each_site{i_ion}(:, 1), ...
            ion_nn_cartesian_each_site{i_ion}(:, 2), 80+40*i_ion, 'red');
        hold on;
        scatter(...
           ion_nnn_cartesian_each_site{i_ion}(:, 1), ...
           ion_nnn_cartesian_each_site{i_ion}(:, 2), 100+40*i_ion, 'g');
    end
end

%%
% Determine left/right orientation (v_ij in Eq. 1 of Kane-Mele 2005)
ion_nnn_left_right_v_ij = cell(n_ion_unit_cell, 1);
for i_site = 1 : n_ion_unit_cell
    n_ion_nn_searching = size(ion_nn_cartesian_each_site{i_site}, 1);
    n_ion_nnn_searching = size(ion_nnn_cartesian_each_site{i_site}, 1);
    ion_nnn_left_right_v_ij{i_site} = zeros(n_ion_nnn_searching, 1);
    i_ion_nnn = 1;
    for i_ion = 1 : n_ion_nn_searching
        searching_nn_ion_cartesian = ...
            ion_nn_cartesian_each_site{i_site}(i_ion, :);
        r_nn_to_nnn = ...
            ion_nnn_cartesian_each_site{i_site} ...
            - searching_nn_ion_cartesian;
        distance_nn_to_nnn = sqrt(sum(r_nn_to_nnn.^2, 2));
        nn_to_nnn_connect_index = abs(distance_nn_to_nnn - nn_distance) < 1e-7;
        r_0_to_nn = ...
            repmat( ...
                searching_nn_ion_cartesian - ion_cartesian_all(i_site, :), ...
                [n_ion_nnn_searching, 1]);
        left_right_index = cross(r_nn_to_nnn, r_0_to_nn, 2);
        left_right_index = sign(left_right_index(:, 3));
        ion_nnn_left_right_v_ij{i_site}(nn_to_nnn_connect_index) = ...
            left_right_index(nn_to_nnn_connect_index);
    end
end

%%
% Generate the on-site and NN hopping terms
hopping_orbit_0 = (1 : 2*n_ion_unit_cell)';
hopping_orbit_R = (1 : 2*n_ion_unit_cell)';
hopping_r = zeros(2*n_ion_unit_cell, 3);
hopping_R = zeros(2*n_ion_unit_cell, 3);
hopping_t = zeros(2*n_ion_unit_cell, 1);
for i_site_0_including_spin = 1 : n_ion_unit_cell*2
    if mod(i_site_0_including_spin, 2) == 1
        hopping_t(i_site_0_including_spin, 1) = hopping_on_site_A;
    else
        hopping_t(i_site_0_including_spin, 1) = hopping_on_site_B;
    end
end

% Add nearest-neighbor (NN) hopping terms
% For graphene, NN hopping only occurs between opposite sublattices (A ↔ B)
for i_site_0_including_spin = 1 : n_ion_unit_cell*2
    i_site_0 = mod(i_site_0_including_spin-1, 2) + 1;
    i_site_R = mod(i_site_0_including_spin, 2) + 1;
    if i_site_0_including_spin == 1
        i_site_R_including_spin = 2;
    elseif i_site_0_including_spin == 2
        i_site_R_including_spin = 1;
    elseif i_site_0_including_spin == 3
        i_site_R_including_spin = 4;
    elseif i_site_0_including_spin == 4
        i_site_R_including_spin = 3;
    end
    n_hopping_nn = size(ion_nn_R_each_site{i_site}, 1);
    hopping_orbit_0 = [ ...
        hopping_orbit_0; 
        i_site_0_including_spin .* ones(n_hopping_nn, 1)];
    hopping_orbit_R = [ ...
        hopping_orbit_R;
        i_site_R_including_spin .* ones(n_hopping_nn, 1)];
    hopping_t = [ ...
        hopping_t;
        hopping_nn_t .* ones(n_hopping_nn, 1)];
    hopping_R = [ ...
        hopping_R;
        ion_nn_R_each_site{i_site_0}];
    hopping_r = [ ...
        hopping_r;
        ion_nn_R_each_site{i_site_0} ...
            + ion_position_all(i_site_R) - ion_position_all(i_site_0)];
end

%%
% Generate next-nearest-neighbor (NNN) hopping terms
for i_site_0_including_spin = 1 : n_ion_unit_cell*2
    spin_index = ((i_site_0_including_spin < n_ion_unit_cell + 1e-5) - 0.5) * 2;
    i_site_R_including_spin = i_site_0_including_spin;
    i_site_0 = mod(i_site_0_including_spin-1, 2) + 1;
    i_site_R = mod(i_site_0_including_spin-1, 2) + 1;
    n_hopping_nnn = size(ion_nnn_R_each_site{i_site}, 1);
    hopping_orbit_0 = [ ...
        hopping_orbit_0; 
        i_site_0_including_spin .* ones(n_hopping_nnn, 1)];
    hopping_orbit_R = [ ...
        hopping_orbit_R;
        i_site_R_including_spin .* ones(n_hopping_nnn, 1)];
    hopping_t = [ ...
        hopping_t;
        (1i .* spin_index .* ion_nnn_left_right_v_ij{i_site_0} ...
            .* hopping_nnn_kane_mele + hopping_nnn_t) ...
            .* ones(n_hopping_nnn, 1)];
    hopping_R = [ ...
        hopping_R;
        ion_nnn_R_each_site{i_site_0}];
    hopping_r = [ ...
        hopping_r;
        ion_nnn_R_each_site{i_site_0} ...
            + ion_position_all(i_site_R) - ion_position_all(i_site_0)];
end

%%
% Save results to file
hopping = table( ...
    hopping_orbit_0, hopping_orbit_R, hopping_t, hopping_R, hopping_r, ...
    'VariableNames', {'orbit_0', 'orbit_R', 't', 'R', 'r'});

% Store hopping parameters and lattice information in a structure
wan_basis = struct();
wan_basis.n_band = 2 * n_ion_unit_cell;
wan_basis.n_hopping = size(hopping, 1);
wan_basis.hopping = hopping;
wan_basis.is_contain_lattice_inf = true;
wan_basis.lattice_a = lattice_a;
wan_basis.lattice_b = lattice_b;

% Save the structure to 'wan_basis.mat'
save wan_basis wan_basis
