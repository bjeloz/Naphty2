function [ dr, r ] = VectorPBC( x1, x2, box_top )

%  Input: x1 and x2 are positions of two atoms or ghostatoms (centers).
%         box_top is the topology of the considered box
% Output: dr corresponds to absolute length of vector r, taking pbc into account

% get box topology data
T = box_top.T;
a = box_top.a;
b = box_top.b;
c = box_top.c;

% calculate distance vector r
r = (x2 - x1);

% r has to be a column vector
sr = size(r);
if sr(1) < sr(2)
    r = r';
end

% r in fractional coordinates
rFrac = T*r;

% pbc for arbitrary box shape
if rFrac(1) > 0.5
    r = r - a;
elseif rFrac(1) < -0.5
    r = r + a;
end
if rFrac(2) > 0.5
    r = r - b;
elseif rFrac(2) < -0.5
    r = r + b;
end
if rFrac(3) > 0.5
    r = r - c;
elseif rFrac(3) < -0.5
    r = r + c;
end

% get the size of distance vector:
dr = sqrt(dot(r,r));

end