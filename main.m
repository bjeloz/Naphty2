clear all;close all;clc;

% System specific data
sys_data.ato_mol = 18;  % number of atoms per molecule

% center of mass calculation
center_indices = [8;9];  % 
%center_indices = 8;  % 

% plot related data
step_size = 0.002;   % histogram stepsize

% read gro file
[gro_ideal_cry, box_dim_ideal_cry] = READ_GRO('data/NAP192Face1m10_T290_NPT_B-equil.gro');  % READ_GRO
%[gro_ideal_cry, box_dim_ideal_cry] = READ_GRO('data/NAPHTA11_5x5x5_supercell.gro');  % READ_GRO

% radial distribution
if false
    
    box_top = FractionalCoordinates(box_dim_ideal_cry);
    
    centers_cry_ideal = Centers(gro_ideal_cry(:,3:5), sys_data, box_top, center_indices); % gro_ideal_cry corresponds to atoms, 

    [X, CDF_cry_ideal] = RadialDistribution(centers_cry_ideal, box_top, step_size);
    
    leg_cry_ideal = plot(X, smooth(CDF_cry_ideal,10), 'Linewidth', 1);
    %leg_cry_ideal = plot(X, CDF_cry_ideal, 'Linewidth', 1);
    grid on
    hold on
    
end


% cumulative radial distribution
if false
    
    box_top = FractionalCoordinates(box_dim_ideal_cry);
    
    centers_cry_ideal = Centers(gro_ideal_cry(:,3:5), sys_data, box_top, center_indices); % gro_ideal_cry corresponds to atoms, 

    [X, CDF_cry_ideal] = CumulativeRadialDistribution(centers_cry_ideal, box_top, step_size);
    
    leg_cry_ideal = plot(X, CDF_cry_ideal, 'Linewidth', 1);
    grid on
    hold on
    
end


% Angular distribution
if false
    
    vector_indices = [8;9];
    
    step_size = 0.02;
    
    box_top = FractionalCoordinates(box_dim_ideal_cry);

    vectors_ideal = Vectors(gro_ideal_cry(:,3:5), sys_data, box_top, vector_indices);
    
    [X, pTheta_cry_ideal] = AngleDistribution(vectors_ideal, step_size);

    leg_cry_ideal = plot(X, smooth(pTheta_cry_ideal, 5), 'Linewidth', 1);
    %leg_cry_ideal = plot(X, pTheta_cry_ideal, 'Linewidth', 1);
    %axis([0 pi/2 0 3])
    grid on
    hold on
    
end