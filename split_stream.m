function [I,Q] = split_stream(mapped)
%copies I column to separate I vector, and Q column to the Q vector
    I=mapped(:,2);
    Q=mapped(:,3);
end