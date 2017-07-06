function [ X, gR ] = RadialDistribution( centers, box_top, step_size )
    % Calculation of the averaged radial distribution function for all N molecules


    a = box_top.a;
    b = box_top.b;
    c = box_top.c;
    V = box_top.V;

    radius = min( [0.5*a(1), 0.5*b(2), 0.5*c(3)] );  % take the minimum of values to calculate the radius

    N = length(centers);

    % Density:
    rho = N / V;

    % histogram of the radial distribution function:
    histogram = zeros(ceil(radius/step_size) + 1, 1);

    for i=1:N-1
        for j=i+1:N

            % get length of intermolecular vector
            [dr,~] = VectorPBC(centers(i,:), centers(j,:), box_top);

            if dr < radius && dr > 0
                % update histogram:
                binIndex = ceil(dr/step_size);
                histogram(binIndex) = histogram(binIndex) + 1;
            end
        end
    end

    gR = histogram(1:length(histogram)-1);
    for i=1:length(gR)
        volShell = 4/3 * pi * ((i*step_size)^3 - ((i-1)*step_size)^3);
        gR(i) = gR(i) / (volShell * rho * N/2);
    end

    % independent variable of histogram, length dr
    X = 0:step_size:(length(gR)-1)*step_size;
    X = X + step_size/2;

end

