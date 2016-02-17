function yy = lagInterPolSyms(x, xp)
    n = length(xp);
    
    for k = 1:n
        temp = 1;
        for i = 1:n
           if(i == k)
               continue;
           end
           temp = temp * ( (x-xp(i)) / (xp(k)-xp(i)));
        end
        yy(k) = {temp};
    end

end