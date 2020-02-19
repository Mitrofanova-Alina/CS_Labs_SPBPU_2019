function [s1, s2, s3, s4] = split_sector_simple(sector, R)
    ind = fix(length(sector) / 2); 

    first_sec = sector(1 : ind);
    second_sec = sector(ind + 1 : length(sector));
    [~ , first_ind] = max(R(first_sec));
    [~ , second_ind] = max(R(second_sec));
   
    second_ind = second_ind + length(first_sec);
    
    len_1 = first_ind;
    len_2 = length(first_sec);
    len_3 = second_ind;
    len_4 = length(sector);
    
    s1 = sector(1 : len_1);    
    s2 = sector(len_1 + 1 : len_2);
    s3 = sector(len_2 + 1 : len_3);
    s4 = sector(len_3 + 1: len_4);
end