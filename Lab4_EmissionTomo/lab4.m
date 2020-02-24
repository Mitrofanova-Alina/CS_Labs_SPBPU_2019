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
%% plot separatrix and magnetic axis
figure()
grid on
hold on
plot(RBDRY, ZBDRY, "ks", 'LineWidth', 1.5);
plot(magnet_axis(1), magnet_axis(2), "r*", 'LineWidth', 1.5);
[x_lim, y_lim] = get_plot_lim(RBDRY, ZBDRY);
xlim(x_lim)
ylim(y_lim)
axis equal
title("separatrix")
legend("separatrix", "magnetic axis")
%%
num_sectors = 8;
num_circle = 6;
cur_cut_y = RBDRY;
cur_cut_z = ZBDRY;
cur_magnet_axis = magnet_axis;

[R_segments_arr, Z_segments_arr, lines_start, lines_end] = get_web_grid(cur_cut_y, cur_cut_z, cur_magnet_axis, num_sectors, num_circle);
elements = create_element_from_grid(R_segments_arr, Z_segments_arr, lines_start, lines_end);
%% plot grid of splitting of separatrix
figure()
grid on
hold on
for i = 1 : length(elements)
    element_.draw(elements(i), "b", i);
end
axis equal
title("grid of splitting separatrix")
legend("grid")
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
detector.horizontal_step = point_.create((detector.end.x - detector.start.x) / 16 , (detector.end.y - detector.start.y) / 16 ); 
%% plot Tokamak view top
figure()
grid on
hold on
axis equal
phi = -pi : pi / 360 : pi;
R = max(RBDRY);
r = min(RBDRY);
plot(R * cos(phi), R * sin(phi), "b", 'LineWidth', 1.5); 
plot(r * cos(phi), r * sin(phi), "r", 'LineWidth', 1.5);

x = detector.start.x;
y = detector.start.y;

for i = 1 : 16
    A = detector_.get_col_pos(detector, i);
    x = A.x : 0.01 : 0.8;
    B = detector.aperture_pos;
    [k, b] = get_line(A, B);
    plot(x, k*x + b, "k--", 'LineWidth', 1)
    text_R = R * 1.2;
    text_x = (-2 * b * k + sqrt(4 * b ^ 2 * k ^ 2 - 4 * (k ^ 2 + 1) * (b ^ 2 - text_R ^ 2))) / (2 *(k ^ 2 + 1));
    text_y = k * text_x + b;
    text(text_x, text_y, num2str(i));
end

xlim([-0.8, 0.8])
ylim([-0.8, 0.8])

title('Tokamak view top')
xlabel('x, m')
ylabel('y, m')
%%
for cut_ind = 1 : 16
    figure()
    grid on
    hold on
    axis equal
    H = detector_.get_plane(detector, cut_ind);
    title(strcat("H = ", num2str(H)));
    element_num = length(elements);
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
    
    for ray_ind = 1 : 16
        [k, b, det_pos, apper_pos] = detector_.get_ray(detector, cut_ind, ray_ind);
        x = det_pos.x : 0.01 : 0.6;
        plot(x, k * x + b, "k")
        for t = 1 : length(cut_elements)
            element_.draw(cut_elements(t), "g")
            [hord, intersection] = element_.get_hord(cut_elements(t), k, b);          
            for g = 1 : length(intersection)
                plot(intersection(g).x, intersection(g).y, "or");
            end
        end
        text_x = 0.6;
        text_y = text_x * k + b;
        text(text_x, text_y, num2str(ray_ind));
    end    
end
%% plot location of detector pixels in its plane
figure()
grid on
hold on
%axis equal
spd_sizes = [0.88, 1.23] * 1e-3;

y_square_up = [0, spd_sizes(1), spd_sizes(1), 0, 0];
y_square_down = [0, 0, spd_sizes(1), spd_sizes(1), 0];
z_square_up = [0, 0, spd_sizes(2), spd_sizes(2), 0];
z_square_down = [0, - spd_sizes(2), - spd_sizes(2), 0, 0];
color = "k";
point_color = ".k";
for j = 0 : 7
    for i = 0 : 7
        x = (spd_sizes(1) + j * sum(spd_xy_step) + y_square_up);
        y = (- spd_z_step - spd_sizes(2)) / 2 + i * (- spd_z_step) + z_square_up;
        plot(x, y, color);
        x = x(1 : 4);
        x = sum(x) / 4;
        y = y(1 : 4);
        y = sum(y) / 4;
        plot(x, y, point_color);
        x = spd_sizes(1) + spd_xy_step(1) + j * sum(spd_xy_step) + y_square_up;
        y = (- spd_z_step - spd_sizes(2)) / 2 + i * (- spd_z_step) + z_square_up;
        plot(x, y, color);
        x = x(1 : 4);
        x = sum(x) / 4;
        y = y(1 : 4);
        y = sum(y) / 4;
        plot(x, y, point_color);
        x = spd_sizes(1) + j * sum(spd_xy_step) + y_square_down;
        y = (spd_z_step + spd_sizes(2)) / 2 + i * spd_z_step + z_square_down;
        plot(x, y, color);
        x = x(1 : 4);
        x = sum(x) / 4;
        y = y(1 : 4);
        y = sum(y) / 4;
        plot(x, y, point_color);
        x = spd_sizes(1) + spd_xy_step(1) + j * sum(spd_xy_step) + y_square_down;
        y = (spd_z_step + spd_sizes(2)) / 2 + i * spd_z_step + z_square_down;
        plot(x, y, color);
        x = x(1 : 4);
        x = sum(x) / 4;
        y = y(1 : 4);
        y = sum(y) / 4;
        plot(x, y, point_color);
    end
end

title('location of detector pixels in its plane');
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
figure
imagesc(hord_matrix);
title('Matrix');
colormap('Bone');
cm = colormap;
colormap(flipud(cm));
colorbar;