function [ cellList, cellNums, cellMols] = createCellList_tric( centers, box_dim, rcut )
% creates the cell list for the speedup of CV calculation
% r_cut corresponds to the list cutting distance and not the CV r_cut one

[T, a, b, c] = fractionalCoordinates(box_dim);

cellNums = [floor(a(1)/rcut), floor(b(2)/rcut), floor(c(3)/rcut)];

meshx = linspace(-0.0001, 1, cellNums(1)+1);
meshy = linspace(-0.0001, 1, cellNums(2)+1);
meshz = linspace(-0.0001, 1, cellNums(3)+1);

%particleList = zeros(length(centers), 4);
cellList = zeros(length(centers), 1);

for i=1:length(centers)
    rFrac = T*centers(i,:)';
    
    % pbc:
    cell_abc = zeros(3,1);
    for j=1:3
        if rFrac(j) < 0
            rFrac(j) = rFrac(j) + 1;
        elseif rFrac(j) > 1
            rFrac(j) = rFrac(j) - 1;
        end
    end
    
    % get a,b,c-inices of cell:
    cell_abc(1) = sum(rFrac(1)>meshx(:));
    cell_abc(2) = sum(rFrac(2)>meshy(:));
    cell_abc(3) = sum(rFrac(3)>meshz(:));
    
%     cellNr = cell_abc(1) + (cell_abc(2)-1)*cellNums(1) ...
%         + (cell_abc(3)-1)*cellNums(1)*cellNums(2);
    
    cellNr = index3Dto1D(cell_abc, cellNums);
    
    %particleList(i,:) = [ centers(i,:) , cellNr ];
    cellList(i) = cellNr;
end

% List of all molecules per cell
sortedList = sort(cellList);
t = 1;
t_old = 0;
maxMolsPerCell = 0;
for i=1:length(cellList)
    if sortedList(i)==t_old
        t = t+1;
    else
        t = 1;
    end
    if t>maxMolsPerCell
        maxMolsPerCell = t;
    end
    t_old = sortedList(i);
end 
cellMols = zeros(cellNums(1)*cellNums(2)*cellNums(3), maxMolsPerCell+1);
for i=1:length(cellList)
    n = sum(cellMols(cellList(i),:)~=0);
    cellMols(cellList(i),n+2) = i;
end
for i=1:length(cellMols)
    cellMols(i,1) = sum(cellMols(i,:)~=0);
end

end

