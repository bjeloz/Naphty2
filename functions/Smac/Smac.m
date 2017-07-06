function [ s ] = Smac( gro_ato, box_dim, sys_data, smac_data )
% Calculation of the Smac collective variable

% CV specifications
center_indices = smac_data.center_indices;
vector_indices = smac_data.vector_indices;
rcut = smac_data.rcut;
ncut = smac_data.ncut;
phi0 = smac_data.phi0;
a = smac_data.a;
b = smac_data.b;

% Rational switching function
rationalSWF = @(x, x0, a, b) (1-(x/x0).^a) ./ (1-(x/x0).^b);

% Gaussian function
bell = @(x, mu, sigma) exp(-(x-mu).^2/(2*sigma^2));

% Vector angle function
vectorAngle = @(v1, v2) acos(dot(v1, v2) ./ (sqrt(dot(v1,v1)).*sqrt(dot(v2,v2))));

% Get simulation box properties
box_top = FractionalCoordinates(box_dim);

% Get predefined centers
centers = Centers( gro_ato, sys_data, box_top, center_indices );
N = length(centers);

% Get triclinic cell list
[cellList, cellNums] = createCellList_tric(centers, box_dim, 1.2*rcut);

s = zeros(N,1);

verletList = zeros(N, 10);
for i=1:N
    neighborList = cellNeighbors(i, centers, cellList, cellNums);

    % get vectors for angles
    v0 = Vectors(gro_ato, sys_data, box_top, vector_indices);
    v_i = Vectors(gro_ato, sys_data, box_top, vector_indices);

    phi_i = zeros(length(neighborList),1);
    for k=1:length(phi0)
        phi_ij = vectorAngle(repmat(v0, length(neighborList), 1)', v_i');
        phi_ij = min(phi_ij, pi-phi_ij);
        phi_i = phi_i + min(1, bell(phi_ij, phi0(k,1), phi0(k,2))');
    end

    f_i = zeros(length(neighborList),1);
    for j=1:length(f_i)
        r_ij = VectorPBC(centers(i,:), centers(neighborList(j),:), box_top);
        f_i(j) = rationalSWF(r_ij, rcut, a, b);
        if f_i(j)>0.5
            verletList(i,1) = verletList(i,1) + 1;
            verletList(i, verletList(i,1)+1) = neighborList(j);
        end
    end

    n_i = sum(f_i);
    rho_i = rationalSWF(n_i, ncut, -a, -b);
    
    s(i) = rho_i/n_i * sum(f_i.*phi_i);
end

end