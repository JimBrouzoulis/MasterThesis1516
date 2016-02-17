function [ ] = el3ddraw( ex,ey,ez )

nel = size(ex,2);
hold on;
for ie = 1:nel
    
   for it = 1:8 
     cel(it,1:3) = [ex(it,ie), ey(it,ie), ez(it,ie)];  
   end
   
   order = [1 2 3 4 1 5 8 7 6 5 8 4 3 7 6 2];
   drawx = cel(order,1);
   drawy = cel(order,2);
   drawz = cel(order,3);
   
   plot3(drawx,drawy,drawz,'k'); 
   
end
hold off;
end

