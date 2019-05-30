function Cell = Matrix2Cell(Matrix)
    num = size(Matrix,3);
    for i =1:num
        Cell{1,i}=Matrix(:,:,i);       
    end
end