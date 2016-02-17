function [edof,coord,ex,ey,ez,dof,nel,ndofs,nno,...
    side1nodes, side2nodes, side3nodes, side4nodes, side5nodes] = cubeMesher(lx,ly,lz,nelx,nely,nelz)
% side1: back, 
% side2: front, 
% side3: left
% side4: right
% side5: top
%                         *
%                      ******
%                  ****   *  **
%               ****      *   **
%            ****         *     **
%         ***             *      **
%       **                *       ***
%       ***               *         **
%       * **              *           **
%       *   **            *          ***
%       *    **           *       ***  *
%       *      **         *    ***     *
%       *       **        *****        *
%       *         **    ***            *
%       *          *****  *            *
%       *            *    *            *
%       *            *    *            *
%       *            *    *            *
%       *            *    *            *
%       *            *    *            *
%       *            *    *            *
%       *            *  *****          *
%       *            ***    **         *
%       *        *****        **       *
%       *     ***    *         **      *
%       *  ***       *           **    *
%       ***          *            **   *
%       **           *              ** *
%         **         *               ***
%          ***       *                **
%            **      *             ***
%             **     *         ****
%               **   *      ****
%                **  *   ****
%                  ******
%                    *

nelno = 8;

elx = lx/(nelx);
ely = ly/nely;
elz = lz/nelz;
xx = 0:elx:lx;
yy = 0:ely:ly;
zz = 0:elz:lz;

nno = (nelx+1)*(nely+1)*(nelz+1);
ndofs = nno*3;
dof = reshape(1:ndofs,3,nno)';
nel = nelx*nely*nelz;

% ellayout = reshape(1:nel, nelx,nely,nelz);
nodelayout = reshape(1:nno, nelx+1,nely+1,nelz+1);

icel = 0;
for iz=0:nelz-1
    for iy=0:nely-1
        for ix = 0:nelx-1
            icel = icel+1;
            temp = nodelayout((1:2) + ix,(1:2) +iy,(1:2) + iz);
            temp = reshape(temp,1,8,1); 
            temp(3:4) = temp([4 3]); temp([7 8]) = temp([8 7]);
            mesh(icel,:) = temp;
            edof(icel,:) = reshape(dof(temp,:)',1,nelno*3);
        end
    end
end
edof = edof';
mesh = mesh';

cn = 0;
for iz=1:nelz+1
    for iy=1:nely+1
        for ix = 1:nelx+1
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





