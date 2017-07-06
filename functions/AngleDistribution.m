function [ X, angleDistribution ] = AngleDistribution( vectors, step_size )

N = length(vectors);

histogram = zeros(ceil(pi/step_size), 1);

for i=1:N-1
    for j=i+1:N
        vi = vectors(i,:);
        vj = vectors(j,:);
                
        angle = acos(dot(vi, vj) / (norm(vi)*norm(vj)));
        
        if angle>pi/2
            angle=pi-angle;
        end
        
        % update histogram:
        binIndex = ceil(angle/step_size);
        histogram(binIndex) = histogram(binIndex) + 1;
    end
end

angleDistribution = histogram ./ (sum(histogram)*step_size);

X = 0:step_size:(length(angleDistribution)-1)*step_size;


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

