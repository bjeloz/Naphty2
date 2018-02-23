function [ s ] = GaussSmac( gro_ato, box_dim, sys_data, smac_data )
    % Calculation of the Smac collective variable

    % CV specifications
    center_indices = smac_data.center_indices;
    vector_indices = smac_data.vector_indices;
    r_cut   = smac_data.r_cut;
    r_cell  = smac_data.r_cell;
    n_cut   = smac_data.n_cut;
    phi_0   = smac_data.phi_0;
    sigma_0 = smac_data.sigma_0;
    r_0     = smac_data.r_0;
    sigma_r = smac_data.sigma_1;
    a       = smac_data.a;
    b       = smac_data.b;

    % Rational switching function
    rationalSWF = @(x, x0, a, b) (1-(x/x0).^a) ./ (1-(x/x0).^b);

    % Gaussian function
    bell = @(x, phi, sigma) exp(-(x-phi).^2/(2*sigma^2));

    % Vector angle function
    vectorAngle = @(v1, v2) acos(dot(v1, v2) ./ (sqrt(dot(v1,v1)).*sqrt(dot(v2,v2))));
    
    % Get simulation box properties
    box_top = FractionalCoordinates(box_dim);

    % Get predefined centers
    centers = Centers( gro_ato, sys_data, box_top, center_indices );
    N = length(centers);

    % Get triclinic cell list
    [cellList, cellNums] = createCellList_tric(centers, box_top, r_cell);

    s = zeros(N,1);

%     verletList = zeros(N, 100);
    for i=1:N
        neighborList = cellNeighbors(i, cellList, cellNums);
        %neighborList = (1:N)';
        
        % get vectors for angles
        v_0 = VectorsCellList(gro_ato, i, sys_data, box_top, vector_indices);
        v_i = VectorsCellList(gro_ato, neighborList, sys_data, box_top, vector_indices);

        Gauss_ij = zeros(length(neighborList),1);
        for j=1:length(Gauss_ij)
            for k=1:length(phi_0)
                phi_ij = vectorAngle(repmat(v_0, length(neighborList(j)), 1)', v_i(j,:)');
                phi_ij = min(phi_ij, pi-phi_ij);
                Gauss_phi_ij = min(1, bell(phi_ij, phi_0(k), sigma_0(k))');
            
                r_ij = VectorPBC(centers(i,:), centers(neighborList(j),:), box_top);
                Gauss_r_ij = min(1, bell(r_ij, r_0(k), sigma_r(k))');
                Gauss_r_ij = Gauss_r_ij*rationalSWF(r_ij, r_cut, a, b);
                Gauss_ij(j) = Gauss_ij(j) + Gauss_r_ij*Gauss_phi_ij;
                if Gauss_ij(j) > 0.2
                    Gauss_ij(j) = 1;
                else
                    Gauss_ij(j) = 0;
                end
            end
        end
        
        s(i) = sum(Gauss_ij);
%         f_i = zeros(length(neighborList),1);
%         for j=1:length(f_i)
%             r_ij = VectorPBC(centers(i,:), centers(neighborList(j),:), box_top);
%             f_i(j) = rationalSWF(r_ij, r_cut, a, b);
%         end

        %n_i = sum(f_i);
        %rho_i = rationalSWF(n_i, n_cut, -a, -b);

        %s(i) = rho_i/n_i * sum(f_i.*phi_i);
        
    end

end