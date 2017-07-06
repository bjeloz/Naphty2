function [ vectors ] = VectorsCellList( gro_ato, list_indices, sys_data, box_top, vector_indices )
% same as Vectors.m, but the algorithm allows for picking certain molecules
% from a given cell of the cell List

    % get number of atoms per molecule
    ato_mol = sys_data.ato_mol;

    % get total number of atoms and molecules
    ato_tot = length(gro_ato);
    mol_tot = ato_tot / ato_mol;

    % if nargin<8
    %     vector_indices = (1:nr_of_mols)';
    % end

    % initialize list of vectors
    vectors = zeros(length(list_indices), 3);
    
    % calculate vectors
    size_vector_indices = length(vector_indices);
    if size_vector_indices < 2
        
        disp('Not enough entries in vector_indices')
        
    elseif size_vector_indices == 2
        
        count=0;
        for ii=0:ato_mol:ato_tot-1
            if any(ii/ato_mol+1 == list_indices)
                count = count+1;
                [~, r] = VectorPBC(gro_ato(ii+vector_indices(1),:), gro_ato(ii+vector_indices(2),:), box_top);
                vectors(count, :) = r;
            end
        end

    elseif size_vector_indices == 3
        
        disp('3 vector_indices are not allowed')
        
    elseif size_vector_indices == 4
        
        count=0;
        for ii=0:ato_mol:ato_tot-1
            if any(ii/ato_mol+1 == list_indices)
                count = count+1;
                [~, r1] = VectorPBC(gro_ato(ii+vector_indices(1),:), gro_ato(ii+vector_indices(2),:), box_top);
                m1 = gro_ato(ii+vector_indices(1),:) + 0.5*r1';
                [~, r2] = VectorPBC(gro_ato(ii+vector_indices(3),:), gro_ato(ii+vector_indices(4),:), box_top);
                m2 = gro_ato(ii+vector_indices(3),:) + 0.5*r2';
                [~, r3] = VectorPBC(m1, m2, box_top);
                vectors(count, :) = r3;
            end
        end

    else
        
        disp('Too many values in vecor_indices')
    
    end
end

