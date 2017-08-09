function [ s1 ] = NNsmooth( s, verletList )
% This function smooths 
% 

s1 = zeros(size(s));
for i=1:length(s)
    n_i = verletList(i,1);
    s1(i) = s(i) / (n_i+1);
    for j=2:n_i+1
        s1(i) = s1(i) + s(verletList(i,j)) / (n_i+1);
    end
end

end

