clear variables
syms x y z real

xx = lagInterPolSyms(x, [-1 1])
yy = lagInterPolSyms(y, [-1 1])
% zz = lagInterPolSyms(z, [-1 -0.5 0 0.5 1])%[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
zz = lagInterPolSyms(z, [-1 -1/3 1/3 1])

itr = 1;
for iz = 1:length(zz)
    for iy = 1:length(yy)
        for ix = 1:length(xx)
            func{itr} = xx{ix}*yy{iy}*zz{iz};
            itr = itr+1;
        end
    end
end

hej(x,y,z) = func{1}
hej(x,y,z) = func
N(x,y,z) = hej; 
%Specific:

B = [diff(N(x,y,z),x); diff(N(x,y,z),y); diff(N(x,y,z),z)]


% for iz = [-1 0 1]
%     for iy = [-1 1]
%         for ix = [-1 1]
%             N(ix,iy,iz)
%         end
%     end
% end
    
    
    
    