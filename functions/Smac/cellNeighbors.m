function [ neighborList ] = cellNeighbors( partIndex, centers, particleList, cellNum )
% Create neighbor list

partCellNr = particleList(partIndex);

[index3D] = index1Dto3D(partCellNr, cellNum);
cellX = index3D(1);
cellY = index3D(2);
cellZ = index3D(3);

% Get cell enums of all 27 cells around and with current cell
neighborCells = zeros(27,1);
count = 1;
for i=-1:1
    for j=-1:1
        for k=-1:1
            neighborCells(count) = mod(cellX-1+i,cellNum(1))+1 ...
                + (mod(cellY-1+j,cellNum(2)))*cellNum(1) ...
                + (mod(cellZ-1+k,cellNum(3)))*cellNum(1)*cellNum(2);
            count = count + 1;
        end
    end
end

% Discard redundancy for small cellNum values (less than 3)
neighborCells = sort(neighborCells);
last = 0;
for i=1:length(neighborCells)
    if neighborCells(i) == last
        neighborCells(i) = 0;
    else
        last = neighborCells(i);
    end
end

% Loop over all particles in particleList and check if they belong to one
% of the cells in neighbourCells
isNeighbor = zeros(length(particleList), 1);
for i=1:length(particleList)
    
    if sum(particleList(i)==neighborCells(:)) == 1 && i~=partIndex
        isNeighbor(i) = 1;
    end
    
end

% Construct a list of all particles contained in the cells of
% neighbourCells
neighborList = zeros(sum(isNeighbor), 1);
count = 1;
for i=1:length(particleList)
    if isNeighbor(i) == 1
        neighborList(count, 1) = i;
        count = count + 1;
    end
end
end

