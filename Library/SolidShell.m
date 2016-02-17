classdef SolidShell < handle
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ir = IntegrationRule;
        
        dispInterp = InterpolatorX2Y2Z2;
        stressInterp;% = InterpolatorX2Y2Z3;
        
        submatrices;
        
        Mhat = -1;
        nnoz = -1; % Number of nodes in z-direction
        
        nDispDofs = -1;
        nStrainDofs = -1;
        nStressDofs = -1;
    end
    
    methods
        function obj = SolidShell(ngpx, ngpy, ngpz, nnoz, Mhat)
           % constructor
           
           obj.ir.setupCubeRule(ngpx, ngpy, ngpz);
           obj.Mhat = Mhat;
           
           if(nnoz == 3)
               obj.stressInterp = InterpolatorX2Y2Z3;
           elseif(nnoz == 4)
               obj.stressInterp = InterpolatorX2Y2Z4;
           else
              error('Interpolator not implemented yet.'); 
           end
           
           obj.nStrainDofs = size(Mhat(0,0,0),2);
           obj.nStressDofs = 6*nnoz*2*2;
           obj.nDispDofs = 3*2*2*2; %*nnoz??
           obj.nnoz = nnoz; %Number of nodes in z-direction
           
        end
        
        function [Kout, fout] = computeLinearizedSystem(obj, ex,ey,ez,eq, D)
            
            Ae = zeros(obj.nStressDofs, obj.nDispDofs  );
            Be = zeros(obj.nStressDofs, obj.nStrainDofs);
            Ce = zeros(obj.nStressDofs, obj.nStressDofs);
            De = zeros(obj.nStrainDofs, obj.nDispDofs);
            Ee = zeros(obj.nStrainDofs, obj.nStrainDofs);
            Fe = zeros(obj.nDispDofs  , obj.nStressDofs);
            fe = zeros(obj.nDispDofs  , 1);
            
            for gp = obj.ir.gps
                
                %Gauss coordinates
                lcoords = gp.local_coords;
                
                %N-vector
                Nvec = obj.dispInterp.eval_N(lcoords);
                
                %Derivatives of n-vector
                [dNdx, detJ] = obj.dispInterp.eval_dNdx(lcoords, ex, ey, ez);
                
                %Get N and B matrix used in FEM
                [N, B] = solid8NandBmatrix(Nvec, dNdx);

                %Enanced part
                M = obj.Mhat(lcoords(1), lcoords(2), lcoords(3));

                %Stress part
                NvecSigma = obj.stressInterp.eval_N(lcoords);
                P = obj.stressInterp.createNmatrix(NvecSigma, 6);
                
                %Integrattion
                dV = detJ * gp.weight; 
                
                Ae = Ae + -P'*D*B * dV;
                Be = Be + -P'*D*M * dV;
                Ce = Ce + P'*P    * dV;
                De = De + M'*D*B  * dV;
                Ee = Ee + M'*D*M  * dV;
                Fe = Fe + B'*P    * dV;
                
                fe = fe + N'*eq   * dV;
            end
            
            obj.submatrices.Ae = Ae;
            obj.submatrices.Be = Be;
            obj.submatrices.Ce = Ce;
            obj.submatrices.De = De;
            obj.submatrices.Ee = Ee;
            obj.submatrices.Fe = Fe;
            
            %Static condensate
            O1 = zeros(obj.nStrainDofs,obj.nStressDofs);
            O2 = zeros(obj.nDispDofs, obj.nDispDofs);
            O3 = zeros(obj.nDispDofs, obj.nStrainDofs);
            
            K = [O2, O3, Fe;...
                 De, Ee, O1;... 
                 Ae, Be, Ce]; 
                    
            f = [fe; zeros(obj.nStrainDofs + obj.nStressDofs,1)];
            
            %We have bc on the stress variables beta, handle it:
            %All dofs
            alldofs = 1:(obj.nDispDofs + obj.nStrainDofs + obj.nStressDofs);
            
            %Bc on simga
            %Bc för utbreddlast på konsolbalk, Txy noll på över och under
            %rand
%             sigmaBc = (obj.nStrainDofs+obj.nDispDofs) + [5:6:24, 53:6:obj.nStressDofs]';
            
            %Bc för Txy noll på över och under-rand, + Szz noll på under
            %rand
            sigmaBc = (obj.nStrainDofs+obj.nDispDofs) + [5:6:24, 53:6:obj.nStressDofs, 3:6:24,  51:6:obj.nStressDofs]';
            sigmaBc = [sigmaBc, sigmaBc*0];
            sigmaBc(end-3:end, 2) = [-40];
            
            nSigmaPredescibed = size(sigmaBc,1);
            
            %Free dofs
            freedofs = setdiff(alldofs,sigmaBc(:,1));
            
            %New matrices
            newK = K(freedofs,freedofs);
            newf = f(freedofs) - K(freedofs, sigmaBc(:,1))*sigmaBc(:,2);     
            
            %Static condenstation
            alldofs = 1:(obj.nDispDofs + obj.nStrainDofs + obj.nStressDofs - nSigmaPredescibed);
            i = 1:obj.nDispDofs;
            d = setdiff(alldofs,i);

            Kout = newK(i,i) - newK(i,d)*inv(newK(d,d))*newK(d,i);
            fout = newf(i) - newK(i,d)*inv(newK(d,d))*newf(d);
            
            %Static condenstation
%             alldofs = 1:(obj.nDispDofs + obj.nStrainDofs + obj.nStressDofs);
%             i = 1:obj.nDispDofs;
%             d = setdiff(alldofs,i);
% 
%             Kout = K(i,i) - K(i,d)*inv(K(d,d))*K(d,i);
%             fout = f(i);
            
        end
        
        function beta = computeStress(obj,ex,ey,ez,a,D)
            De = obj.submatrices.De;
            Ee = obj.submatrices.Ee;
            
            alpha = -Ee\(De*a);
            
            Ae = obj.submatrices.Ae;
            Be = obj.submatrices.Be;
            Ce = obj.submatrices.Ce;
            
            %----
            allSigmaDofs = 1:obj.nStressDofs;
%             sigmaLockedDofs = [5:6:24, 53:6:obj.nStressDofs]';
            sigmaLockedDofs = [5:6:24, 53:6:obj.nStressDofs, 3:6:24,  51:6:obj.nStressDofs]';
%             sigmaBc = [sigmaBc, sigmaBc*0];
            sigmaBc = [sigmaLockedDofs, sigmaLockedDofs*0];
            sigmaBc(end-3:end, 2) = [-40];
            
            sigmaFreeDofs = setdiff(allSigmaDofs,sigmaLockedDofs);
            
            beta = zeros(obj.nStressDofs,1);
            beta(sigmaLockedDofs) = sigmaBc(:,2);
            
            tempF = -(Ae*a + Be*alpha);
            beta(sigmaFreeDofs) = Ce(sigmaFreeDofs, sigmaFreeDofs)\(tempF(sigmaFreeDofs) - Ce(sigmaFreeDofs,sigmaLockedDofs)*sigmaBc(:,2));
        end
        
        function R = computeR(obj,a,ex,ey,ez,D) % residual
            error('Not ready for use, Update the Interpolator-changes');
            
            %Init
            A = zeros(obj.nStrainDofs,obj.nDispDofs);
            C = zeros(obj.nStrainDofs);
            
            for gp = obj.ir.gps
                
                %Get gauss coordinates
                lcoords = gp.local_coords;
                
                %Get N-vector
                Nvec = obj.interp.eval_N(lcoords, 2);
                
                %Get B-vectorr
                [dNdx, detJ] = obj.interp.eval_dNdx(lcoords, ex, ey, ez, 2);
                
                %Get N and B matrix for Fem
                [~, B] = solid8NandBmatrix(Nvec, dNdx);
                
                %Enanced part
                M = obj.Mhat(lcoords(1), lcoords(2), lcoords(3));

                %Integration
                dV = detJ * gp.weight;
                A = A + M'*D*B * dV;
                C = C + M'*D*M * dV;                
                %fext = fext + N'*eq * dV; % body forces
            end
            
            % solve for alpha: A'*a + C*alpha = 0
            alpha = -C\(A*a);
            
            % Determine stresses from a least square fit
            %  int( P'*sig ) - int(P'*P)*beta = 0
            L    = zeros(obj.nStressDofs, obj.nStressDofs);
            fsig = zeros(obj.nStressDofs, 1);
            Ae = zeros(obj.nStressDofs, obj.nDispDofs  );
            Be = zeros(obj.nStressDofs, obj.nStrainDofs);
            Fe = zeros(obj.nDispDofs  , obj.nStressDofs);
            
            for gp = obj.ir.gps
                %Gauss coordinates
                lcoords = gp.local_coords;
                
                %N-vecor
                Nvec = obj.interp.eval_N(lcoords, 2);
                
                %Stress matrix
                NvecSigma = obj.interp.eval_N(lcoords, 3);
                P = obj.interp.createNmatrix(NvecSigma, 6);
                
                %Enhanced matrix
                M = obj.Mhat(lcoords(1), lcoords(2), lcoords(3));
                
                %Derivatives of N-vec
                [dNdx, detJ] = obj.interp.eval_dNdx(lcoords, ex, ey, ez, 2);
                
                %N and B matrix for FEM
                [~, B] = solid8NandBmatrix(Nvec, dNdx);
                
                %Compatible and enhances part
                epsC = B*a;
                epsA = M*alpha;
                
                %True stress, so to speak
                sig = D * (epsC + epsA);
                
                %Integrate
                dV = detJ * gp.weight;
                L = L + P'*P * dV;
                fsig = fsig + P'*sig * dV;
                
                %Matrices for post-processing stresses
                Ae = Ae - P'*D*B * dV;
                Be = Be - P'*D*M * dV;
                Fe = Fe + B'*P   * dV;
            end
            
            %Silly bc to test if the problem is solvable
             allSigmaDofs = 1:obj.nStressDofs;
            sigmaLockedDofs = [5:6:24, 53:6:obj.nStressDofs]';
            sigmaFreeDofs = setdiff(allSigmaDofs,sigmaLockedDofs);
            
            beta = solveq(L,fsig, [sigmaLockedDofs, sigmaLockedDofs*0]);%,stress_bc);
            
            % compute residual
            R = zeros(obj.nDispDofs,1);
            R_ref = R;
            for gp = obj.ir.gps
                %Gauss coordinates
                lcoords = gp.local_coords;
                
                %N-vector
                Nvec = obj.interp.eval_N(lcoords, 2);
                
                %Derivatives of N-vector
                [dNdx, detJ] = obj.interp.eval_dNdx(lcoords, ex, ey, ez, 2);
                
                %Fem B and N matrices
                [~, B] = solid8NandBmatrix(Nvec, dNdx);
                
                %ENhanced part
                M = obj.Mhat(lcoords(1), lcoords(2), lcoords(3));
                
                %Compatible and assumed strains
                epsC = B*a;
                epsA = M*alpha;
                
                %For reference
                sig_old = D * (epsC + epsA);
                
                %Matrix for stress
                NvecSigma = obj.interp.eval_N(lcoords, 3);
                P = obj.interp.createNmatrix(NvecSigma, 6);
                sig = P*beta;
                
                %A simple check
                err = abs(sig-sig_old)./(abs(sig_old)+1e-12);
                if err > 1e-3
                  error('Note: LSF of stresses differ from compatible ones') 
                end
                
                %Integration
                dV = detJ * gp.weight;
                R = R + B' * sig * dV;
%                 R = R + B' * sig_old * dV;
            end
            
            obj.submatrices.Ae = Ae;
            obj.submatrices.Be = Be;
            obj.submatrices.Ce = L;
            obj.submatrices.De = A;
            obj.submatrices.Ee = C;
            obj.submatrices.Fe = Fe;
            
            
%                 abs(R-R_ref)
%             % pure EAS
%             R = zeros(obj.nDispDofs,1);
%             for gp = obj.ir.gps
%                 
%                 lcoords = gp.local_coords;
%                 Nvec = obj.interp.eval_N(lcoords);
%                 [dNdx, detJ] = obj.interp.eval_dNdx(lcoords, ex, ey, ez);
%                 [N, B] = solid8NandBmatrix(Nvec, dNdx);
%                 
%                 
%                 M = Mhat(lcoords(1), lcoords(2), lcoords(3));
%                 epsC = B*a;
%                 epsA = M*alpha;
%                 sig = D * (epsC + epsA);
%    
%                 %Integration
%                 dV = detJ * gp.weight; % delta V
%                 R = R + B'*sig * dV;
%             
%             end
            
            
 
        end
       
        function Knum = computeNumTangent(obj, R0, ae, ex, ey, ez, D, delta)
            % delta - size of numerical perturbation
            neldofs = length(R0);
            Knum = zeros(neldofs);
            
            for i=1:neldofs
                apert = ae;
                apert(i) = apert(i) + delta;
                Rpert = obj.computeR(apert, ex, ey, ez, D);
                Knum(:,i) = (Rpert-R0) / delta; 
            end
        end
    end
    
end

