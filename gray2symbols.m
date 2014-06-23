function [ outvec ] = gray2symbols( invec, map )
% mapping I, Q pairs to decimal values


for i = 1:length(invec(:,1))
    for j = 1:length(map)
        if invec(i,:)==map(j,2:3);
        outvec(i)=map(j);
    end

end
end
end

