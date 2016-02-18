
%Post processing of "consol-beam"

%Max displacement
maxDisp = max(abs(a));

%Area moment of inertia
Iy = (ly*lz^3/12); 

%Analytical values from Euler-bernoulli
tauxz = @(x,z,b,h,L,q) 3/2/(b*h)*(1 - (z/(h/2)).^2)*q;   %q = [N]
sigzz = @(x,z,b,h,L,q) zeros(size(z));
eb_maxdisp = @(P,ly, lx,E,Iy) abs(P*lx^3/3/E/Iy); % P = [N)

%Compare deflections between FEM and EB
eb_maxdisp = eb_maxdisp(P,ly,lx,E,Iy);
fprintf('EulerBernoulli: %.10f, SolidElement: %.10f\n',eb_maxdisp,maxDisp);


%STRESSES
for iel = 1:nel
    
    %Element Stresses
    elStress(:,iel) = el(iel).computeStress(ex(:,iel)',ey(:,iel)',ez(:,iel)', ed(:,iel), D);
    
    %Calculate the global coordinates from the local
    plotZZlocal = linspace(-1, 1, 10);
    for iz = 1:length(plotZZlocal)
        temp = el(iel).dispInterp.eval_globalCoords([-1,-1, plotZZlocal(iz)], ex(:,iel)',ey(:,iel)',ez(:,iel)');
        plotZZglobal(iz) = temp(3);
    end
    plotXXlocal  = [-1 0 1];
    for ix = 1:length(plotXXlocal)
        temp = el(iel).dispInterp.eval_globalCoords([plotXXlocal(ix),-1,-1], ex(:,iel)',ey(:,iel)',ez(:,iel)');
        plotXXglobal(ix) = temp(1);
    end
    
    
    %Stresses
    for ix = 1:length(plotXXlocal)
        
        %Calucalete stresses
        for iz = 1:length(plotZZlocal)
            
            N = el(iel).getStressPmatrix([plotXXlocal(ix) 0 plotZZlocal(iz)]);
            
            tempStress = N*elStress(:,iel);
            
            %NOTE: Here you change what stress to plot
            %
            plotStress(iz) = tempStress(5); %<---------------
        end
        
        %Analtycial values
        tauxzAnalytical = tauxz(plotXXglobal(ix),plotZZglobal - lz/2,ly,lz,lx,P); %PPP **** LYYYYYY
        sigzzAnalytical = sigzz(plotXXglobal(ix),plotZZglobal - lz/2,ly,lz,lx,P); %PPP **** LYYYYYY
        
        %Plot
        tf = subplot(1,3,ix);
        plot(plotStress,plotZZglobal);  hold on;
        plot(tauxzAnalytical,plotZZglobal,'r')
%         plot(sigzzAnalytical,plotZZglobal,'r')
        
        title(sprintf('Element: %i, \n Local x-coord: [%0.3f]',iel, plotXXlocal(ix)));
        

        legend('FEM','Euler-Bernoulli')
    end
    
    %Break here, so you can step through the stresses in different elements
    keyboard;
    clf
end