function [result] = HeurMinCond(A_inf, A_sup, iter_number)
    % find size of matrix
    m = size(A_inf,1); 
    n = size(A_inf,2);   
    
    % setting the number of random throws in the implemented algorithm
    if(nargin >= 3)
        NN = iter_number;
    else
        NN = 10;
    end
    
    % initializing corner matrices for A
    Matr1 = ones(m,n); 
    Matr2 = ones(m,n);
    
    % initialize MinCond- the minimum of condition numbers of point matrices 
    % contained in the specified interval matrix A
    MinCond = Inf; 
    for k = 1:NN 
        % randomly generate an integer matrix EPM from zeros and
        % of units, the same size as A ( random interval
        % of integer values is specified by the first randi argument)
        EPM = randi([0,1],m,n); 
          
        % generate angular matrices that are diagonally opposite
        % to each other, according to the "diagonal" heuristic
        % optimization method
        for i = 1:m
            for j = 1:n
                if EPM(i,j) == 0 
                    Matr1(i,j) = A_inf(i,j); 
                    Matr2(i,j) = A_sup(i,j); 
                else 
                    Matr1(i,j) = A_sup(i,j);
                    Matr2(i,j) = A_inf(i,j);                 
                end 
            end
        end 
        
        % find the condition number of the obtained
        % of angular matrices , adjusting the minimum estimate
        c1 = cond(Matr1,2); 
        c2 = cond(Matr2,2); 
        if MinCond > c1 
            MinCond = c1; 
        end
        if MinCond > c2 
            MinCond = c2; 
        end
        
    end
    
    result = MinCond;
end