function [ out ] = FractionalCoordinates( box_dim )
    % calculates:
    % transformation matrix T
    % box side vectors a, b, c in Cartesian coordinates
    % and the box volume V
    %
    % T*[x,y,z] provides fractional coordinates [xFrac,yFrac,zFrac] for point [x,y,z]
    %
    % fractional coordinates are necessary for the calcuation of pbc
    % check www.cchem.berkeley.edu/chem195/radial_distribution_8m.html

    a = [ box_dim(1); box_dim(4); box_dim(5) ];
    b = [ box_dim(6); box_dim(2); box_dim(7) ];
    c = [ box_dim(8); box_dim(9); box_dim(3) ];

    % Calculate box angles:
    da = sqrt(sum(dot(a,a)));
    db = sqrt(sum(dot(b,b)));
    dc = sqrt(sum(dot(c,c)));
    cos_alpha = dot(b,c) / (db*dc);
    cos_beta = dot(a,c) / (da*dc);
    cos_gamma = dot(a,b) / (da*db);
    sin_gamma = sqrt(1 - cos_gamma^2);

    % Triclinic box volume:
    V = da*db*dc*sqrt(1-cos_alpha^2-cos_beta^2-cos_gamma^2 + ...
        2*cos_alpha*cos_beta*cos_gamma);

    % Transformation matrix T; [xFrac;yFrac;zFrac] = T * [x;y;z]
    T = [1/da, -cos_gamma/(da*sin_gamma), ...
        db*dc*(cos_alpha*cos_gamma-cos_beta)/(V*sin_gamma); ...
        0, 1/(db*sin_gamma), da*dc*(cos_beta*cos_gamma-cos_alpha)/(V*sin_gamma); ...
        0, 0, da*db*sin_gamma/V];

    out.T = T;
    out.a = a;
    out.b = b;
    out.c = c;
    out.V = V;
end

