function [ tauxz, sigzz, eb_maxdisp] = analyticalStuff( input_args )

%Konsolbalk med punktlast.
% tauxz = @(x,z,b,h,L,q) 3/2/(b*h)*(1 - (z/(h/2)).^2)*q;   %q = [N]
% sigzz = @(x,z,b,h,L,q) 0;
% eb_maxdisp = @(P,ly, lx,E,Iy) abs(P*lx^3/3/E/Iy); % P = [N)

%Konsolbalk med utbredd last.
tauxz = @(x,z,b,h,L,q) 3/2/(b*h)*(1 - (z/(h/2)).^2)*q*(L-x); %q ? [N/m]
sigzz = @(x,z,b,h,L,q) (q*(h-z).*(h + 2*z).^2)/(2*(b*h)*h^2);
eb_maxdisp = @(P,ly,lx,E,Iy) abs((P*ly)*lx^4)/(8*E*Iy); %P = [N/m]

end

