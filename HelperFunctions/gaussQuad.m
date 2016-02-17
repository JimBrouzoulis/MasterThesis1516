function [points,weights]=gaussQuad(nx2,ny2,nz2)

ny = 0; nz = 0;
ngp =1;

if(nargin < 1)
    error('Not enough input arguments.');
else
    [p1,w1] = gausQuadTable(nx2);
    ngp = ngp*nx2;
    nx = nx2;
    if(nargin > 1)
        [p2,w2] = gausQuadTable(ny2);
        ngp = ngp*ny2;
        ny = ny2;
        if(nargin > 2)
            [p3,w3] = gausQuadTable(nz2);
            ngp = ngp*nz2;
            nz = nz2;
            if(nargin > 3)
                error('maybe in the future');
            end
        end
    end
end

dim = nargin;
points = zeros(dim,ngp);
weights = zeros(1,ngp);
itr = 1;
for ii = 1:nx
    ptemp(1) = p1(ii);
    w1temp = w1(ii);
    for jj = 1:ny
        ptemp(2) = p2(jj);
        w2temp = w2(jj);
        for kk = 1:nz
            ptemp(3) = p3(kk);
            w3temp = w3(kk);
            points(:,itr) = ptemp;
            weights(:,itr) = w1temp*w2temp*w3temp;
            itr = itr+1;
        end
        
        if(nz == 0)
            points(:,itr) = ptemp;
            weights(:,itr) = w1temp*w2temp;
            itr = itr+1;
        end
    end
    
    if(ny == 0)
        points(:,itr) = ptemp;
        weights(:,itr) = w1temp;
        itr = itr+1;
    end
end

end

