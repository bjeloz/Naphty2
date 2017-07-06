function [ centers ] = Centers( gro_ato, sys_data, box_top, center_indices )
    % gives as output the center of each naphthalene molecule 
    % the input atoms come from a .gro file

    % get number of atoms per molecule
    ato_mol = sys_data.ato_mol;
    
    % get total number of atoms and molecules
    ato_tot = length(gro_ato);
    mol_tot = ato_tot / ato_mol;
    
    % initialize data for calculation of center values
    centers = zeros(mol_tot, 3);
    count = 0;

    % calculate centers
    size_center_indices = length(center_indices);  % center_enum gives the atom positions from which the center is calculated
    if size_center_indices == 1
        
        for i = center_indices(1):ato_mol:ato_tot
            count = count + 1;
            centers(count,:) = gro_ato(i,:);
        end
        
    elseif size_center_indices == 2
        
        delta_center_enum = center_indices(2)-center_indices(1);
        for i = center_indices(1):ato_mol:ato_tot
            count = count + 1;
            [~, r] = VectorPBC(gro_ato(i,:), gro_ato(i+delta_center_enum,:), box_top);
            centers(count,:) = gro_ato(i,:)+0.5*r';
        end
        
    else
        
        disp('Too many values in center_indices')
        
    end
end

