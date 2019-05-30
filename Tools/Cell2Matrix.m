function [Cov_Matrix,Log_Matrix]= Cell2Matrix(Cell)
    num = length(Cell);
    for i =1:num
        Cov_Matrix(:,:,i)=Cell{1,i};  
        Log_Matrix(:,:,i)=logm(Cell{1,i});
    end
end