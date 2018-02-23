function [ X, cosAngleDistribution ] = CosAngleDistribution( vectors, step_size )

% Function AngleDistribution needs following functions
% for the calculation of the input:
% -> FractionalCoordinates.m
% -> Vectors.m
% -> VectorPBC.m

N = length(vectors);

histogram = zeros(ceil(2/step_size), 1);

for i=1:N-1
    for j=i+1:N
        vi = vectors(i,:);
        vj = vectors(j,:);
                
        cosangle = dot(vi, vj) / (norm(vi)*norm(vj));
        
        if cosangle>0
            cosangle=-cosangle;
        end
        
        cosangle = cosangle + 1;
        
        % update histogram:
        binIndex = round(cosangle/step_size + 0.5);
        histogram(binIndex) = histogram(binIndex) + 1;
    end
end

cosAngleDistribution = histogram ./ (sum(histogram)*step_size);

X = -1:step_size:(length(cosAngleDistribution)/2-1)*step_size;


% trial to obtain pdf(cos(phi)) vs. cos(phi) distribution
% N = length(vectors);
% 
% histogram = zeros(ceil(2/step_size), 1); % cos(phi) spans from -1 to 1
% size_histogram = length(histogram);
% 
% for i=1:N-1
%     for j=i+1:N
%         vi = vectors(i,:);
%         vj = vectors(j,:);
%                 
%         cosangle = dot(vi, vj) / (norm(vi)*norm(vj));
%         
%         if cosangle > 0
%             cosangle = cosangle-1;
%         end
%         
%         % update histogram:
%         binIndex = ceil(cosangle/step_size+size_histogram/2);
%         histogram(binIndex) = histogram(binIndex) + 1;
%     end
% end
% 
% angleDistribution = histogram ./ (sum(histogram)*step_size);
% 
% X = -1:step_size:1-step_size;

end

