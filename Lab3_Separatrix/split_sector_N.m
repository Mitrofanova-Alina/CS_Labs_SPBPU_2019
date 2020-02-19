function [result] = split_sector_N(sector, N)
    result = {};
    cur_step = ceil(length(sector) / N);
    cur_left = 1;
    cur_right = cur_step;
    reduse = ceil(length(sector));
    for i = 1 : N
        if(reduse == 0)
            return
        end
        tmp = sector(cur_left : cur_right);
        result = {result{:} tmp};
        reduse = reduse - length(tmp);
        cur_left = cur_right + 1;
        cur_step = ceil(reduse / (N - i));
        cur_right = cur_right + cur_step;  
    end 
end