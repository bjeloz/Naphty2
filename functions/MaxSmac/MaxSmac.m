function [ s ] = maxsmac( atoms, box_dim, rcut, ncut )
% calculates maxsmac CV, angles are compared for their correlation
% atoms corresponds to all atom coordinates of the *.xtc or *.gro files

%fermi = @(x, m) (exp(x-m)+1).^(-1);
rationalSWF = @(x, x0, a, b) (1-(x/x0).^a) ./ (1-(x/x0).^b);
bell = @(x, mu, sigma) exp(-(x-mu).^2/(2*sigma^2));
vectorAngle = @(v1, v2) acos(dot(v1, v2) ...
    ./(sqrt(dot(v1,v1)).*sqrt(dot(v2,v2))));

a = 15;
b = 30;

[T, av, bv, cv] = fractionalCoordinates(box_dim);   % 
centers = centers_8_9(atoms, box_dim, T, av, bv, cv);
N = length(centers);
apm = length(atoms)/N;
[cellList, cellNums] = createCellList_tric(centers, box_dim, rcut); % create cell list

s = zeros(N,1);  % initialize size of CV of each molecule

verletList = zeros(N, 10);  % internal verlet list
for i=1:N
    neighborList = cellNeighbors(i, centers, cellList, cellNums); % look for each molecule its neighbours in the cell list
    
    % get vectors for angles
    v0 = vectors_1_2_5_7(atoms, i, apm, box_dim, T, av, bv, cv);
    v_i = vectors_1_2_5_7(atoms, neighborList, apm, box_dim, T, av, bv, cv);
%     v02 = vectors_8_9(atoms, i, apm, box_dim, T, av, bv, cv);
%     v_i2 = vectors_8_9(atoms, neighborList, apm, box_dim, T, av, bv, cv);

%     n1 = cross(repmat(v01, length(neighborList), 1)', v_i1');
%     n2 = cross(repmat(v02, length(neighborList), 1)', v_i2');
    phi_i = real(vectorAngle(repmat(v0, length(neighborList), 1)', v_i')');
    % Mirror angle>90? at 90?:
    phi_i = min(phi_i, pi-phi_i);
    
    f_i = zeros(length(neighborList),1);
    % f_i = zeros(length(centers),1);
    for k=1:length(f_i)     % creates verlet list
        r_ij = normPBC(centers(i,:), centers(neighborList(k),:), box_dim, T, av, bv, cv); % calculates distance to the neighboring molecules
    %     r_ij = normPBC(centers(mol_nr,:), centers(j,:), box_dim);
        f_i(k) = rationalSWF(r_ij, rcut, a, b);
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
    rho_i = rationalSWF(n_i, ncut, -a, -b);
    
%     if rho_i/n_i*PHI_i>1
%         rho_i/n_i*PHI_i
%     end
    %s(i) = -1/100*log(sum(exp(-100*[rho_i/n_i*PHI_i, 0.9999])));%* sum(f_i.*PHI_i));
    s(i) = rho_i/n_i*PHI_i;
end

s = NNsmoth(s, verletList);     % CV vector for each molecule

end

