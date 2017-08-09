%clear all;close all;clc;

% System specific data
sys_data.ato_mol = 18;  % number of atoms per molecule

% center of mass calculation
%center_indices = [8;9];  % 
center_indices = 8;  % 

% plot related data
step_size = 0.002;   % histogram stepsize

% read gro file
%[gro_ideal_cry, box_dim_ideal_cry] = READ_GRO('data/angletests/NAP320ETH713Face00m1_T300_NPT.r0_0.0ps.gro');  % READ_GRO
%[gro_ideal_cry, box_dim_ideal_cry] = READ_GRO('data/angletests/NAP320ETH713Face00m1_T300_NPT.r0_1000.0ps.gro');  % READ_GRO

%[gro_ideal_cry, box_dim_ideal_cry] = READ_GRO('data/NAP192Face1m10_T290_NPT_B-equil.gro');  % READ_GRO
[gro_ideal_cry, box_dim_ideal_cry] = READ_GRO('data/NAP473-NPT-T300-init.gro');  % READ_GRO
%[gro_ideal_cry, box_dim_ideal_cry] = READ_GRO('data/NAPHTA11_5x5x5_supercell.gro');  % READ_GRO

% radial distribution
if false
    
    box_top = FractionalCoordinates(box_dim_ideal_cry);
    
    centers_cry_ideal = Centers(gro_ideal_cry(:,3:5), sys_data, box_top, center_indices); % gro_ideal_cry corresponds to atoms, 

    [X, CDF_cry_ideal] = RadialDistribution(centers_cry_ideal, box_top, step_size);
    
    %leg_cry_ideal = plot(X, smooth(CDF_cry_ideal,10), 'Linewidth', 1);
    leg_cry_ideal = plot(X, CDF_cry_ideal, 'Linewidth', 1);
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
    %vector_indices = [1;2;5;7];
    
    step_size = 0.02;
    
    box_top = FractionalCoordinates(box_dim_ideal_cry);

    vectors_ideal = Vectors(gro_ideal_cry(:,3:5), sys_data, box_top, vector_indices);
    
    [X, pTheta_cry_ideal] = AngleDistribution(vectors_ideal, step_size);

    %leg_cry_ideal = plot(X, smooth(pTheta_cry_ideal, 5), 'Linewidth', 1);
    leg_cry_ideal = plot(X, pTheta_cry_ideal, 'Linewidth', 1);
    %axis([0 pi/2 0 3])
    grid on
    hold on
    
end

% Smac
if true
    
    % Vector [8;9]
    smac_data.center_indices = [8];
    smac_data.vector_indices = [8; 9];
    smac_data.r_cut = 0.65;
    smac_data.r_cell = 1.2*smac_data.r_cut;  % cut-off distance of cell list
    smac_data.n_cut = 4;
    %smac_data.phi_0 = [0; 0.99];  % we just need to write the vector angle values below pi/2
    smac_data.phi_0 = [0; 1.00];  % we just need to write the vector angle values below pi/2
    smac_data.sigma_0 = [0.25; 0.20];
    smac_data.a = 12;
    smac_data.b = 24;

%     % Vector [1;2;5;7]
%     smac_data.center_indices = [8; 9];
%     smac_data.vector_indices = [1; 2; 5; 7];
%     smac_data.r_cut = 0.65;
%     smac_data.r_cell = 1.2*smac_data.r_cut;  % cut-off distance of cell list
%     smac_data.n_cut = 5;
%     smac_data.phi_0 = [0; 0.43];  % we just need to write the vector angle values below pi/2
%     smac_data.sigma_0 = [0.15; 0.15];
%     smac_data.a = 15;
%     smac_data.b = 30;

    step_size = 0.02;

    gsmac_ideal = Smac(gro_ideal_cry(:,3:5), box_dim_ideal_cry, sys_data, smac_data);
    histogram = zeros(ceil(1.01/step_size),1);
    for i=1:length(gsmac_ideal)
        binIndex = ceil(gsmac_ideal(i)/step_size+0.0001);
        histogram(binIndex) = histogram(binIndex) + 1.00;
    end
    gsmac_ideal_dist = histogram/(sum(histogram)*step_size);
    X = linspace(0, 1.01, length(histogram));
    %leg_cry_ideal = plot(X, smooth(gsmac_ideal_dist,5), 'Linewidth', 1);
    leg_cry_ideal = plot(X, gsmac_ideal_dist, 'Linewidth', 1);
    hold on
    
end