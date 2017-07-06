function [ X, gR_int ] = CumulativeRadialDistribution( centers, box_top, step_size )
    % Calculation of the averaged cumulative radial distribution function for all n molecules

    %T = box_top.T;
    a = box_top.a;
    b = box_top.b;
    c = box_top.c;

    radius = min( [0.5*a(1), 0.5*b(2), 0.5*c(3)] );  % take the minimum of the values to calculate the radius

    N = length(centers);
    
    % histogram of cumulative radial distribution function:
    histogram = zeros(ceil(radius/step_size) + 1, 1);  % ceil rounds content to the nearest integer greater than or equal to that content

    for i=1:N-1
        for j=i+1:N

            % get length of intermolecular vector:
            [dr,~] = VectorPBC(centers(i,:), centers(j,:), box_top);

            if dr < radius && dr > 0 % check if molecule pair is in the defined region of the CRDF
                % update histogram:
                binIndex = ceil(dr/step_size);
                histogram(binIndex) = histogram(binIndex) + 1;
            end
        end
    end
    
    gR_int = histogram(1:length(histogram)-1);  % do not take last value of histogram into account (why?)
    gR_int = gR_int/(N/2);  % divide by half of the number of considered molecules 
                            % (since we are using all possible molecule pairs for the calculation of gR_int)

    for i=2:length(gR_int)
        gR_int(i) = gR_int(i-1) + gR_int(i);
    end
    
    % independent variable of histogram, length dr
    X = 0:step_size:(length(gR_int)-1)*step_size;

end

