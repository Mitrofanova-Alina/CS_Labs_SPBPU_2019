%%
clear all
x_lim_left = -0.015;
x_lim_right = 0.825;  % empirically selected constants for axis scale alignment
%% read gfile
[flux,RBDRY,ZBDRY,NBDRY,R,Z,time,rdim,zdim]=gfile_extractor_1t(34363,000162,65);
%%
[arr, ind_arr] = min(flux);
[flux_min, min_j] = min(arr);
min_i = ind_arr(min_j);
%% plot separatrix
figure()
grid on
hold on
plot([RBDRY, RBDRY(1)], [ZBDRY, ZBDRY(1)], "ks", 'LineWidth', 1.5);
magnet_axis = [R(min_j), Z(min_i)];
plot(magnet_axis(1), magnet_axis(2), "r*", 'LineWidth', 1.5);
xlim([x_lim_left x_lim_right]) 
title("separatrix")
legend("separatrix", "magnetic axis")
%% calculate radius of curvature separately for the first and last point, taking into account the closure of the curve
a = [RBDRY(NBDRY), ZBDRY(NBDRY)];
b = [RBDRY(1), ZBDRY(1)];
c = [RBDRY(2), ZBDRY(2)];      
R = find_curv_radius(a, b, c);

for i = 2 : NBDRY - 1
    a = [RBDRY(i - 1),ZBDRY(i - 1)];
    b = [RBDRY(i),ZBDRY(i)];
    c = [RBDRY(i + 1),ZBDRY(i + 1)];
    tmp = find_curv_radius(a, b, c);
    R = [R, tmp];
end

a = [RBDRY(NBDRY - 1),ZBDRY(NBDRY - 1)];
b = [RBDRY(NBDRY),ZBDRY(NBDRY)];
c = [RBDRY(1),ZBDRY(1)];
tmp = find_curv_radius(a, b, c);
R = [R, tmp];
%% plot radius of curvature
figure()
grid on
hold on
plot(R, "b-", 'LineWidth', 1.5);
xlabel('radius');
ylabel('angle');
title("radius of curvature")
%% calculate smoothed radius of curivature
R_sm = medfilt1(R);
%% plot radius of curvature
figure()
grid on
hold on
plot(R_sm, "b-", 'LineWidth', 1.5);
xlabel('radius');
ylabel('angle');
title("smoоthed radius of curvature")
%% calculate grid (gossamer)
N = 5;
init_sector = 1 : NBDRY;
[s1, s2, s3, s4] = split_sector_simple(init_sector, R_sm);
sectors = {s1, s2, s3, s4};
result = {};
for i = 1 : length(sectors)
    tmp_sec = sectors{i};
    tmp = split_sector_N(tmp_sec, N);
    result = {result{:} tmp{:}};
end
sectors = result;

seg = [];
middle = magnet_axis;
R_seg = [];
Z_seg = [];

for i = 1 : length(sectors)
    tmp_sec = sectors{i};
    ind = tmp_sec(1);
    point = [RBDRY(ind), ZBDRY(ind)];
    mid_point = find_center(point, middle);
    R_seg = [R_seg, mid_point(1)];
    Z_seg = [Z_seg, mid_point(2)];
end
%% plot grid
figure()
grid on
hold on
plot([RBDRY, RBDRY(1)], [ZBDRY, ZBDRY(1)], "ks", 'LineWidth', 1.5);
plot(magnet_axis(1), magnet_axis(2), "*r", 'LineWidth', 2);
xlim([x_lim_left x_lim_right]) 
title("grid")


for i = 1 : length(sectors)
    tmp_sector = sectors{i};
    ind = tmp_sector(1);
    point = [RBDRY(ind),ZBDRY(ind)];
    plot([point(1), middle(1)], [point(2), middle(2)], "b-", 'LineWidth', 1.5);
end
for i = 1 : length(R_seg) - 1
    plot([R_seg(i), R_seg(i + 1)], [Z_seg(i), Z_seg(i + 1)], "b-", 'LineWidth', 1.5);
end
plot([R_seg(length(R_seg)), R_seg(1)], [Z_seg(length(R_seg)), Z_seg(1)], "b-", 'LineWidth', 1.5);
plot(magnet_axis(1), magnet_axis(2), "*r", 'LineWidth', 2);