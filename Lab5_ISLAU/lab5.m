%%
clear all
point_ = get_point_interface();
element_ = get_element_interface();
detector_ = get_detector_interface();
%%
[flux, RBDRY, ZBDRY, NBDRY, R, Z, time, rdim, zdim] = gfile_extractor_1t(34363, 000162, 65);
%%
% change the start of the crawl and change the detour of the separatrix to counter-hourly
% since the last point is equal to the first, let's take this into account
%a new beginning - the 32nd element
[RBDRY,ZBDRY] = circle_spin_and_reverse(RBDRY,ZBDRY,NBDRY, 32);

[arr, ind_arr] = min(flux);
[flux_min, min_j] = min (arr);
min_i = ind_arr(min_j);
magnet_axis = [R(min_j), Z(min_i)];
%%
num_sectors = 8;
num_circle = 6;
cur_cut_y = RBDRY;
cur_cut_z = ZBDRY;
cur_magnet_axis = magnet_axis;

[R_segments_arr, Z_segments_arr, lines_start, lines_end] = get_web_grid(cur_cut_y, cur_cut_z, cur_magnet_axis, num_sectors, num_circle);
elements = create_element_from_grid(R_segments_arr, Z_segments_arr, lines_start, lines_end);
%% calculate partition grids of the separatrix in the cross section plane x = H
% H - distance from the center of the tokamak to the section plane
min_x_sep = min(RBDRY);
max_x_sep = max(RBDRY);
H = 0;
h = length(elements);
cut_elements = [];
left_cut_elem = [];
right_cut_elem = [];
for i = 1 : h
    % the result of the intersection of the sector rotation shape and the x = H plane
    [elem, count] = element_.get_cut(elements(i), H);
    if(count == 2)
        left_cut_elem = [left_cut_elem, elem(1)];
        right_cut_elem = [right_cut_elem, elem(2)];
    else
        left_cut_elem = [left_cut_elem, elem]; 
    end
end

cut_elements = [left_cut_elem, right_cut_elem];
%% magic constants
% angle between the pinhole camera direction and the center direction (between 8 and 9 beams)
ang = acos((708^2 + 720^2 - 31^2) / (2 * 708 * 720));
% position of the detector edge (1st column)
spd_start = [0, -0.708];
% position of the 16th column
spd_end = [0.72 * sin(ang), 0.72 * -cos(ang)];
% direction vector of the pinhole camera in the Equatorial plane
spd_vect = (spd_end - spd_start) / norm(spd_end - spd_start);
% step between columns in the detector plane, 2 numbers
spd_xy_step = [2.3375 - 0.88 , 3.81 - 2.3375 + 0.88 ] * 1e-03;
% the center of the detector
pp = spd_start + spd_vect * ((spd_xy_step(1) + spd_xy_step(2)) * 8 + 0.52 * 1e-03) / 2;
% aperture offset from the center of the detector
aperture_xy_offset = 0.0395;
% coordinate of the aperture
aperture_xy = [pp(1) - spd_vect(2) * aperture_xy_offset, pp(2) + spd_vect(1) * aperture_xy_offset];
spd_z_start = (27.52 - 0.49) / 2 * 1e-03;
spd_z_step = -1.72 * 1e-03;
spd_xy = spd_start + spd_vect * (spd_xy_step(2) / 2 + 0.26 * 1e-03);
%% detectors parametres
detector = detector_.create();
detector.start = point_.create(spd_start(1), spd_start(2));
detector.end = point_.create(spd_end(1), spd_end(2));
detector.step = spd_xy_step;
detector.direction = point_.create(spd_vect(1), spd_vect(2));
detector.center = point_.create(pp(1), pp(2));
detector.aperture_offset = 0.0395;
detector.aperture_pos = point_.create(aperture_xy(1), aperture_xy(2)); 
detector.z_start = spd_z_start;
detector.z_step = spd_z_step;
%%
element_num = length(elements);
hord_matrix = zeros(256, element_num);
cur_row = 1;
for cut_ind = 1 : 16
    H = detector_.get_plane(detector, cut_ind);
    cut_elements = [];
    for i = 1 : element_num
        [elem, count] = element_.get_cut(elements(i), H);
        if(H < min_x_sep)
            if(count == 2)
                cut_elements = [cut_elements, elem(2)];
            else
                cut_elements = [cut_elements, elem];
            end
        else
            cut_elements = [cut_elements, elem];
        end
    end
    N = length(cut_elements);
    for ray_ind=1 : 16
        [k, b, det_pos, apper_pos] = detector_.get_ray(detector, 16, ray_ind);
        for t = 1 : N
            [hord, intersection] = element_.get_hord(cut_elements(t), k, b); 
            hord_matrix(cur_row, cut_elements(t).index) = hord_matrix(cur_row, cut_elements(t).index) + hord;
        end
        cur_row = cur_row + 1;
    end
end

figure('Name', 'HORD MATRIX')
hold on
grid on
imagesc(hord_matrix);
title('37000 Hord Matrix');
colormap('Bone');
cm = colormap;
colormap(flipud(cm));
colorbar;
print('-dpng', '-r300', strcat('37000_hord_matrix', '.png')), pwd;
%%
input_file_name = strcat(num2str(37000), "_SPD16x16.mat");
input_data = load(input_file_name);
sign_bb = input_data.sign_bb(:, :, :);
cnt_meas = size(sign_bb, 3);
tp = cell2mat(input_data.Data(1, 2)) * 1e-3;
tz = cell2mat(input_data.Data(2, 2));
t_start = tz;
t_end = t_start + (cnt_meas - 1) * tp;
t_i = t_start : tp : t_end;
%% overview of integrated luminosity
t_cons_start = 125;
t_cons_end = 200;

dt_cons = 1;
start_efit_time_i = t_cons_start : dt_cons : t_cons_end;

B = [];
for start_efit_time = t_cons_start : t_cons_end
    ind = find(abs(t_i - start_efit_time) < tp / 2);
    b = [];
    for i = 16 : -1 : 1
        b = [b; sign_bb(16: -1 : 1, i, ind(1))];
    end
    b = double(b);
    Bnew = sum(b(:));
    B = [B, Bnew];
end
%% plot overview of integrated luminosity
figure('Name', 'LUMINOSITY OVERVIEW')
hold on
grid on
plot([t_cons_start : t_cons_end], B, 'LineWidth', 1);
title(['37000 SPD16x16.mat', ' Sum b']);
xlabel('start efit time');
print('-dpng', '-r300', strcat('37000_over_lum', '.png')), pwd;
%% select the "window" by which we calculate the boundaries of b
b_time_window = 1;
input_time_period = 000162;
b_data = [];
for start_efit_time = input_time_period - b_time_window : input_time_period + b_time_window
    ind = find(abs(t_i - start_efit_time) < tp / 2);
    b = [];
    for i = 16 : -1 : 1
        b = [b; sign_bb(16 : -1 : 1, i, ind(1))];
    end
    b = double(b);
    b_data = [b_data, b];
end
N = length(b_data);
inf_b = zeros(N, 1);
sup_b = zeros(N, 1);
for i = 1 : N
    inf_b(i) = min(b_data(i, :));
    sup_b(i) = max(b_data(i, :));
end
%%
b = (sup_b + inf_b) / 2;
A = hord_matrix;
%% simple decision
simple_dec = inv(A' * A) * A' * b;
lambda = eig(A' * A);
disp(strcat("cond(A) = ", num2str(cond(A))));
disp(strcat("cond(A'A) = ", num2str(cond(A' * A))));
disp(strcat("number of lambda > 0.2: ", num2str(length(find(lambda > 0.2)))));

figure('Name', "EIG A'A")
hold on
grid on
histogram(lambda)
xlabel("lambda")
ylabel("number")
title("37000 Eigenvalue A'A")
print('-dpng', '-r300', strcat('37000_lambda', '.png')), pwd;
%% tolsolvty first
b_sup = sup_b;
b_inf = inf_b;
A_inf = A;
A_sup = A;
[tolmax, argmax, envs, ccode] = tolsolvty(A_inf, A_sup, b_inf, b_sup);
disp(strcat("tolmax = ", num2str(tolmax)));

figure('Name', "FIRST DECISION tolsolvty")
hold on
grid on
x_axis = 1 : length(b_sup);
y_axis = A_sup * argmax;
plot(x_axis, y_axis');
plot(x_axis, b_inf');
plot(x_axis, b_sup');
xlabel("i")
ylabel("value")
title("37000 First decision tolsolvty")
legend("Ax", "inf b", "sup b")
print('-dpng', '-r300', strcat('37000_first_tolsolvty', '.png')), pwd;
%% tolsolvty second
delta_b = - tolmax;
disp(strcat("delta b = ", num2str(delta_b)));
b_sup = b_sup + delta_b;
b_inf = b_inf - delta_b;
[tolmax, argmax, envs, ccode] = tolsolvty(A_inf, A_sup, b_inf, b_sup);
disp(strcat("tolmax = ", num2str(tolmax)));

figure('Name', "SECOND DECISION tolsolvty")
hold on
grid on
x_axis = 1 : length(b_sup);
y_axis = A_sup * argmax;
plot(x_axis, y_axis');
plot(x_axis, b_inf');
plot(x_axis, b_sup');
xlabel("i")
ylabel("value")
title("37000 Second decision tolsolvty (extended b)")
legend("Ax", "inf b", "sup b")
print('-dpng', '-r300', strcat('37000_second_tolsolvty', '.png')), pwd;
%% plot decision
figure('Name', "OBTAINED SOLUTION tolsolvty")
hold on
grid on
plot(argmax, "*", 'LineWidth', 1)
xlabel("i")
ylabel("x_i")
title("37000 Obtained solution by tolsolvty (extended b)")
print('-dpng', '-r300', strcat('37000_obt_sol', '.png')), pwd;
%% plot hist decision
figure('Name', "HISTOGRAM of OBTAINED SOLUTION tolsolvty")
hold on
grid on
histogram(argmax)
xlabel("x_i")
ylabel("number")
title(["37000 Histogram of obtained solution", " by tolsolvty (extended b)"])
print('-dpng', '-r300', strcat('37000_hist_obt_sol', '.png')), pwd;
%% fix radius, different number of iteration
disp("fix radius, different number of iteration")
matrix_radius = 0.1;
A_inf_1 = A * (1 - matrix_radius);
A_sup_1 = A * (1 + matrix_radius);
b_inf_1 = inf_b;
b_sup_1 = sup_b; 
for i = 10 : 10 : 100
    cond_A = HeurMinCond(A_inf_1, A_sup_1, i);
    disp(strcat("rad = ", num2str(matrix_radius), " HeurMinCond(A, ", num2str(i), ") = ", num2str(cond_A)));
end
%% fix number of iteration, different radius
disp("fix number of iteration, different radius")
for rad = 0.1 : 0.05 : 0.5
    A_inf_1 = A * (1 - rad);
    A_sup_1 = A * (1 + rad);
    b_inf_1 = inf_b;
    b_sup_1 = sup_b;
    cond_A = HeurMinCond(A_inf_1, A_sup_1, i);
    disp(strcat("rad = ", num2str(rad), " HeurMinCond(A, ", num2str(i), ") = ", num2str(cond_A)));
end
%% calculate IVE
A_inf_1 = A * 0.9;
A_sup_1 = A * 1.1;
b_inf_1 = inf_b;
b_sup_1 = sup_b;

cond_A = HeurMinCond(A_inf_1, A_sup_1);
[tolmax_1, argmax_1, envs_1, ccode_1] = tolsolvty(A_inf_1, A_sup_1, b_inf_1, b_sup_1);
if (tolmax_1 < 0)
    delta_b_1 = - tolmax_1;
    b_inf_1 = b_inf_1 - delta_b_1;
    b_sup_1 = b_sup_1 + delta_b_1;
    [tolmax_1, argmax_1, envs_1, ccode_1] = tolsolvty(A_inf_1, A_sup_1, b_inf_1, b_sup_1);
end

A_IVE_1 = IVE(A_inf_1, A_sup_1, b_inf_1, b_sup_1, tolmax_1, argmax_1, length(argmax_1));


