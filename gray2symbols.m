function [ outvec ] = gray2symbols( invec, map )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



for i = 1:length(invec)
    for j = 1:length(map)
        if invec(i,:)==map(j,2:3);
        outvec(i)=map(j);
    end

end
end
end

