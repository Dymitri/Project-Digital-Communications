function [ outvec ] = map2gray( invec, map )
%maps vector of symbols to the array outvec containing symbol, I , Q
%according to delivered map arrays

tmpvec = invec*(2.^(size(invec,2)-1:-1:0))';

for i = 1:length(tmpvec)
    for j = 1:length(map)
        if tmpvec(i)==map(j);
        outvec(i, :)=map(j, :);
    end

end
end
end