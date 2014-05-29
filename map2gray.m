function [ outvec ] = map2gray( invec, map )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tmpvec = invec*(2.^(size(invec,2)-1:-1:0))';
%while tmpvec(i)~=map(i)
%    i=i+1;
%    outvec=map(i,:,:);
%end


for i = 1:length(tmpvec)
    for j = 1:length(map)
        if tmpvec(i)==map(j);
        outvec(i, :)=map(j, :);
    end

end
end
end