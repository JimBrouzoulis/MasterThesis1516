function [edof,coord,ex,ey,ez,dof,nel,ndofs,nno,...
    side1nodes, side2nodes, side3nodes, side4nodes, side5nodes] = cubeMesherHigherOrder(lx,ly,lz,nelx,nely,nelz,nnoxel,nnoyel,nnozel,ndofsno)
%Only works proporly for 2x2 nodes in-plane
% side1: back, 
% side2: front, 
% side3: left
% side4: right
% side5: top
%                         7
%                      ******
%                  ****   *  **
%               ****      *   **
%            ****         *     **
%         ***             *      **
%       8*                *       ***
%       ***               *         **
%       * **              *           *6
%       *   **            *          ***
%       *    **           *       ***  *
%       *      **         *    ***     *
%       *       **        *****        *
%       *         **    ***            *
%       *          **5**  *            *
%       *            *    *            *
%       *            *    *            *
%       *            *    *            *
%       *            *    *            *
%       *            *    *            *
%       *            *    3            *
%       *            *  *****          *
%       *            ***    **         *
%       *        *****        **       *
%       *     ***    *         **      *
%       *  ***       *           **    *
%       ***          *            **   *
%       4*           *              ** *
%         **         *               ***
%          ***       *                *2
%            **      *             ***
%             **     *         ****
%               **   *      ****
%                **  *   ****
%                  ******
%                    1

if(nargin < 7)
    ndofsno = 3;
end

nelno = nnoxel*nnoyel*nnozel;%Number of element nodes

elx = lx/(nelx);
ely = ly/nely;
elz = lz/nelz;
xx = 0:elx:lx;
yy = 0:ely:ly;
zz = 0:elz:lz;

totnnox = (nelx*(nnoxel-1)+1);
totnnoy = (nely*(nnoyel-1)+1);
totnnoz = (nelz*(nnozel-1)+1);


xx = linspace(0,lx,totnnox);%0:elx:lx;
yy = linspace(0,ly,totnnoy);
zz = linspace(0,lz,totnnoz);


nno = totnnox*totnnoy*totnnoz;

ndofs = nno*ndofsno;
dof = reshape(1:ndofs,ndofsno,nno)';
nel = nelx*nely*nelz;

% ellayout = reshape(1:nel, nelx,nely,nelz);
nodelayout = reshape(1:nno, totnnox, totnnoy, totnnoz);

icel = 0;
for iz=0:(nnozel-1):totnnoz-nnozel
    for iy=0:(nnoyel-1):totnnoy-nnoyel
        for ix = 0:(nnoyel-1):totnnox-nnoyel
            icel = icel+1;
            temp = nodelayout((1:nnoxel) + ix,(1:nnoyel) +iy,(1:nnozel) + iz);
            temp = reshape(temp,1,nelno,1); 
%             temp(3:4) = temp([4 3]); temp([7 8]) = temp([8 7]);
%             for itt=1:nnozel
%                 temp((3:4) + 4*(itt-1)) = temp(fliplr((3:4) + 4*(itt-1)))
%             end
            mesh(icel,:) = temp;
            edof(icel,:) = reshape(dof(temp,:)',1,nelno*ndofsno);
            
        end
    end
end
edof = edof';
mesh = mesh';

cn = 0;
for iz=1:totnnoz
    for iy=1:totnnoy
        for ix = 1:totnnox
            cn = cn +1;
            coord(cn,:) = [xx(ix), yy(iy), zz(iz)];
        end
    end
end
coord = coord';

xcoord = coord(1,:);
ycoord = coord(2,:);
zcoord = coord(3,:);

ex = xcoord(mesh);
ey = ycoord(mesh);
ez = zcoord(mesh);

if size(ex,1) == 1 & size(ex,2) > 1
   ex = ex'; ey = ey'; ez = ez'; 
end

side1nodes = nodelayout(1,:,:);
side1nodes = side1nodes(:);
side2nodes = nodelayout(end,:,:);
side2nodes = side2nodes(:);

side3nodes = nodelayout(:,1,:);
side3nodes = side3nodes(:);

side4nodes = nodelayout(:,end,:);
side4nodes = side4nodes(:);

side5nodes = nodelayout(:,:,end);
side5nodes = side5nodes(:);
