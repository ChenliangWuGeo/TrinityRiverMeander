function accuS  = streamLineDistance(ds)
    [~,a2] = size(ds);
    for m = 1:a2
        accuS(m) = sum(ds(1:m));
    end
end