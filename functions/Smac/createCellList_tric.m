function [ cellList, cellNums ] = createCellList_tric( centers, box_top, rcell )
    % creates the cell list for the speedup of CV calculation
    % r_cut corresponds to the list cutting distance and not the CV r_cut one

    % get box properties
    T = box_top.T;
    a = box_top.a;
    b = box_top.b;
    c = box_top.c;

    % get number for cells for all three spacial directions
    cellNums = [floor(a(1)/rcell), floor(b(2)/rcell), floor(c(3)/rcell)];  % floor rounds to integer towards negative infinity

    % get mesh of cells in fractional coordinates with equidistant mesh nodes
    % for each dimension:
    meshx = linspace(-0.0001, 1, cellNums(1)+1);
    meshy = linspace(-0.0001, 1, cellNums(2)+1);
    meshz = linspace(-0.0001, 1, cellNums(3)+1);

    %particleList = zeros(length(centers), 4);
    cellList = zeros(length(centers), 1);

    for i=1:length(centers)
        centersFrac = T*centers(i,:)';  % convert center positions from cartesian coordinates to fractional coordinates (for PBC)

        % PBC:
        cell_abc = zeros(3,1);
        for j=1:3
            if centersFrac(j) < 0
                centersFrac(j) = centersFrac(j) + 1;
            elseif centersFrac(j) > 1
                centersFrac(j) = centersFrac(j) - 1;
            end
        end

        % get a,b,c-inices of cell:
        cell_abc(1) = sum(centersFrac(1)>meshx(:));  % sum(value>vector(:)) gives the number of vector elements smaller than given value
        cell_abc(2) = sum(centersFrac(2)>meshy(:));
        cell_abc(3) = sum(centersFrac(3)>meshz(:));

        cellNr = index3Dto1D(cell_abc, cellNums);   % provides 1D vector which numbers the cells of the cell list

        cellList(i) = cellNr;  % 1D cell list

    end

end

