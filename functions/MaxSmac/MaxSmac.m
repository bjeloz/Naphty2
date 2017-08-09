function [ s ] = MaxSmac( gro_ato, box_dim, sys_data, maxsmac_data )
    % calculates maxsmac CV, angles are compared for their correlation
    % atoms corresponds to all atom coordinates of the *.xtc or *.gro files

    % System specifications
    ato_mol = sys_data.ato_mol;

    % CV specifications
    center_indices = maxsmac_data.center_indices;
    vector_indices = maxsmac_data.vector_indices;
    r_cut = maxsmac_data.r_cut;
    r_cell = maxsmac_data.r_cell;
    n_cut = maxsmac_data.n_cut;
    a = maxsmac_data.a;
    b = maxsmac_data.b;

    % Rational switching function
    rationalSWF = @(x, x0, a, b) (1-(x/x0).^a) ./ (1-(x/x0).^b);

    % Gaussian function
    bell = @(x, mu, sigma) exp(-(x-mu).^2/(2*sigma^2));

    % Vector angle function
    vectorAngle = @(v1, v2) acos(dot(v1, v2) ...
        ./(sqrt(dot(v1,v1)).*sqrt(dot(v2,v2))));

    % Get simulation box properties
    box_top = FractionalCoordinates(box_dim);

    % Get predefined centers
    centers = Centers( gro_ato, sys_data, box_top, center_indices );
    N = length(centers);

    % Get triclinic cell list
    [cellList, cellNums] = createCellList_tric(centers, box_top, r_cell); % create cell list

    s = zeros(N,1);  % initialize size of CV of each molecule

    verletList = zeros(N, 10);  % internal verlet list
    for i=1:N
        neighborList = cellNeighbors(i, cellList, cellNums); % look for each molecule its neighbours in the cell list

        % get vectors for angles
        v_0 = VectorsCellList(gro_ato, i, sys_data, box_top, vector_indices);
        v_i = VectorsCellList(gro_ato, neighborList, sys_data, box_top, vector_indices);

        phi_i = real(vectorAngle(repmat(v_0, length(neighborList), 1)', v_i')');
        phi_i = min(phi_i, pi-phi_i); % Mirror angle>90deg at 90deg:

        f_i = zeros(length(neighborList),1);
        % f_i = zeros(length(centers),1);
        for k=1:length(f_i)     % creates verlet list
            r_ij = VectorPBC(centers(i,:), centers(neighborList(k),:), box_top); % calculates distance to the neighboring molecules
            f_i(k) = rationalSWF(r_ij, r_cut, a, b);
            if f_i(k)>0.5       % store nearest neighbors for the maxsmac smooting
                verletList(i,1) = verletList(i,1) + 1;
                verletList(i, verletList(i,1)+1) = neighborList(k);
            end
        end

    %     PHI_i = 0;%zeros(length(neighborList),1);
    %     for k=1:length(phi_i)-1
    %         for l = k+1:length(phi_i)
    %             PHI_i = PHI_i + real(bell(phi_i(k), phi_i(l), 0.175)'*f_i(k)*f_i(l));
    %         end
    %         
    %     end

        PHI_i = 0;%zeros(length(neighborList),1);
        for k=1:length(phi_i)
            PHI_ik = zeros(length(phi_i),1);
            for l = 1:length(phi_i)
                if l~=k
                    %PHI_ik = PHI_ik + exp(300*real(bell(phi_i(k), phi_i(l), 0.05)'*f_i(k)*f_i(l)));
                    PHI_ik(l) = real(bell(phi_i(k), phi_i(l), 0.05)'*f_i(k)*f_i(l));
                end
            end
            %PHI_i = PHI_i + sum(PHI_ik);
            %PHI_i = PHI_i + max(PHI_ik);
            %PHI_i = PHI_i + norm(PHI_ik, 8);
            PHI_i = PHI_i + 1/300*log(sum(exp(PHI_ik*300)));
            %PHI_i = PHI_i + 1/300*log(PHI_ik);
        end

        n_i = sum(f_i);
        %n_i = (n_i)*(n_i-1);
        rho_i = rationalSWF(n_i, n_cut, -a, -b);

    %     if rho_i/n_i*PHI_i>1
    %         rho_i/n_i*PHI_i
    %     end
        %s(i) = -1/100*log(sum(exp(-100*[rho_i/n_i*PHI_i, 0.9999])));%* sum(f_i.*PHI_i));
        s(i) = rho_i/n_i*PHI_i;
    end

    s = NNsmooth(s, verletList);     % CV vector for each molecule

end

