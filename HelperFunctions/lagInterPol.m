function y = lagInterPol(x, xp)
    n = length(xp);
    y = zeros(size(xp));
    
    for k = 1:n
        temp = 1;
        for i = 1:n
           if(i == k)
               continue;
           end
           temp = temp * ( (x-xp(i)) / (xp(k)-xp(i)));
        end
        y(k) = temp;
%         yy(k) = {temp}
    end

end