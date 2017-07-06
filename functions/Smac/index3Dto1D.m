function [ index1D ] = index3Dto1D( index3D, cellNums )
%Index conversion when a 3D array is converted into a 1D array
%   @index3D:   3D index = [ d1, d2, d3 ]
%               domain(d1)=[1,D1], domain(d2)=[1,D2], domain(d3)=[1,D3]
%   @cellNums:  size of 3D array (array of length 3)
%               ==> [ D1, D2, D3 ]
%   @index1D:   1D index
%               domain(index1D)=[1,D1*D2*D3]

index1D = index3D(1) + (index3D(2)-1)*cellNums(1) ...
        + (index3D(3)-1)*cellNums(1)*cellNums(2);

end

