function split_stream(x)
    lx = (length(x));
    half = ceil(lx/2);
    s1 = x(1:half)
    s2 = x(half + 1:lx)
end