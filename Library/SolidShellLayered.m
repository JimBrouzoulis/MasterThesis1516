classdef SolidShellLayered < handle
    %SolidShellLayered Sold-shell element for a laminate with different ply
    %properties.
    
    properties
        
        %Layered integration rule
        lir;% = IntegrationRule;
        
        %Interolationrules
        dispInterp;
        stressInterp;
        stressInterpOrder; %interpolation order for each stresss
        
        submatrices;
        
        %Interpolation matrix for strains
        Mhat = -1;
        
        %Number of dofs
        nDispDofs = -1;
        nStrainDofs = -1;
        nStressDofs = -1;
        nLamStrainDofs = -1;
        
        %Element coords
        ez; ey; ex;
        
        %Number of laminas
        nLam;
        
        %The deterimatant of each layer (always concstant for each element)
        lay_detJ;
    end
    
    methods
        function obj = SolidShellLayered(ngpx, ngpy, ngpz, ex, ey, ez, stressInterp, Mhat, elprop)
                
            %Number of layers
            obj.nLam = elprop.nLam;
            
            %Create integration rules
            obj.lir = LayeredIntegrationRule(elprop.nLam, ngpx, ngpy, ngpz);
            
            %Create Strain interpolation matrices
            obj.Mhat = Mhat;
            
            %Define interpolation for each stress component
            for ii=1:length(stressInterp)
                val = stressInterp(ii);
                if(val == 2)
                    obj.stressInterp{ii} = InterpolatorX2Y2Z2;
                elseif(val == 3)
                    obj.stressInterp{ii} = InterpolatorX2Y2Z3;
                elseif(val == 4)
                    obj.stressInterp{ii} = InterpolatorX2Y2Z4;
                else
                    error('Interpolator not implemented yet.');
                end
            end
            obj.stressInterpOrder = stressInterp;
            
            %Define number of dofs
            obj.nLamStrainDofs = size(Mhat(0,0,0),2);
            obj.nStrainDofs = obj.nLamStrainDofs * obj.nLam;
            obj.nStressDofs = sum(stressInterp*2*2);
            obj.nDispDofs = 3*2*2*2; %*nnoz??
            
            %Save element coords
            obj.ex = ex; obj.ey = ey; obj.ez = ez;
            
            %Interpolation for displacement always the same
            obj.dispInterp = InterpolatorX2Y2Z2;
            
        end
        
        function [Kout, fout] = computeLinearizedSystem(obj,eq,eTrac, elprop)
            
            Ae = zeros(obj.nLamStrainDofs, obj.nDispDofs , obj.nLam);
            Le = zeros(obj.nStressDofs, obj.nLamStrainDofs , obj.nLam);
            Ce = zeros(obj.nLamStrainDofs, obj.nLamStrainDofs , obj.nLam);
            Re = zeros(obj.nStressDofs, obj.nDispDofs);
            Qe = zeros(obj.nStressDofs, obj.nStressDofs);
            Ke = zeros(obj.nDispDofs  , obj.nStressDofs);
            
            fe = zeros(obj.nDispDofs  , 1);
            
            V = 0;
            for ilay = 1:obj.nLam
                
                layD = elprop.Dmatrices(:,:,ilay);
                
                for gp = obj.lir.irs(ilay).gps

                    %Gauss coordinates in local-layer-system and
                    %local-element-system (only need to map the z-coord
                    %since the x and y-coords should be the same 
                    lcoords = gp.local_coords; 
                    ecoords = lcoords;
                    ecoords(3) = obj.lir.getElementGaussCoordinate(ilay,...
                        lcoords(3), elprop.int_coordsL(ilay:ilay+1));
                    
                    %N-vector, eval in element-system
                    Nvec = obj.dispInterp.eval_N(ecoords);

                    %Derivatives of n-vector
                    [dNdx, detJEl] = obj.dispInterp.eval_dNdx(ecoords, obj.ex, obj.ey, obj.ez);

                    % Get detJ for the layer
                    layZ = [[1,1,1,1]*elprop.int_coordsG(ilay), [1,1,1,1]*elprop.int_coordsG(ilay+1)];
                    [~, detJ] = obj.dispInterp.eval_dNdx(ecoords, obj.ex, obj.ey, layZ);
                    
                    % Get N and B matrix used in FEM
                    [N, B] = solid8NandBmatrix(Nvec, dNdx);

                    % Enanced part, eval in layered-system
                    Mi = obj.Mhat(lcoords(1), lcoords(2), lcoords(3));

                    % Stress part, eval in element-system
                    P = [];
                    % Varför inte är inte alla 6 spänningskomponenter med? /JB
                    for is = 1:length(obj.stressInterp)
                        NvecSigma = obj.stressInterp{is}.eval_N(ecoords);
                        % TODO: Mycket långsamt med blkdiag! /JB
                        P = blkdiag(P,NvecSigma);
                    end

                    % Integration
                    dV = detJ * gp.weight;
                    V = V + 1*dV;
                    Ae(:,:,ilay) = Ae(:,:,ilay) + Mi'*layD*B  * dV;
                    Ce(:,:,ilay) = Ce(:,:,ilay) + Mi'*layD*Mi * dV;
                    Le(:,:,ilay) = Le(:,:,ilay) + P'*layD*Mi  * dV;
                    Ke =           Ke           + B'*P        * dV;
                    Re =           Re           + P'*layD*B   * dV;
                    Qe =           Qe           + P'*P        * dV;
                                       
                    fe = fe + N'*eq   * dV;
                end                
            end
            AreaTemp = 4;
            Ntrac = obj.dispInterp.eval_N([0 0 1]);
            Ntrac = obj.dispInterp.createNmatrix(Ntrac, 3);
            
            % Förstår ej 1000-faktorn /JB
            ftrac = Ntrac'* eTrac * AreaTemp* (detJEl*1000); %multiply by a 1000 because detJ is for volume, not for area.
            
            fe = fe + ftrac;
            
            AA = []; LL = []; CC = [];
            for il = 1:obj.nLam
               AA = [AA; Ae(:,:,il)];
               CC = blkdiag(CC, Ce(:,:,il));
               LL = [LL, Le(:,:,il)];
            end
            
            O1 = zeros(obj.nDispDofs, obj.nDispDofs + obj.nStrainDofs);
            O2 = zeros(obj.nStrainDofs, obj.nStressDofs);
            K = [O1, Ke; AA, CC O2; Re, LL, -Qe];
            
            f = [fe; zeros(obj.nStrainDofs + obj.nStressDofs,1)];
%             keyboard;
            %We have bc on the stress variables beta, handle it:
            %All dofs
            alldofs = 1:(obj.nDispDofs + obj.nStrainDofs + obj.nStressDofs);
            
            %%Bc on simga, %InterPolation z=3
            %Bc för utbreddlast på konsolbalk
            %Txy noll på över och under-rand
            [TxybDofs, TxytDofs]  = obj.getStressComponentTopAndBottonDofs(5);
            TxyDofs = [TxybDofs, TxytDofs];
            sigmaBc = (obj.nStrainDofs+obj.nDispDofs) + TxyDofs';
            
            sigmaBc = [sigmaBc, sigmaBc*0];
            nSigmaPredescibed = size(sigmaBc,1);
            
            %Free dofs
            freedofs = setdiff(alldofs,sigmaBc(:,1));
            
            %New matrices, with the known Beta-variables removed
            newK = K(freedofs,freedofs);
            newf = f(freedofs) - K(freedofs, sigmaBc(:,1))*sigmaBc(:,2);
            
            %%%%%% INGA randvilkor:
            %%%
            %
%             newK = K; newf = f;
%             nSigmaPredescibed = 0;
            %
            %%%
            %%%%%
            
            %Static condenstation
            nDofs = obj.nDispDofs + obj.nStrainDofs + obj.nStressDofs;
            alldofs = 1:(nDofs - nSigmaPredescibed);
            i = 1:obj.nDispDofs;
            d = setdiff(alldofs,i);
            
            Kout = newK(i,i) - newK(i,d)*inv(newK(d,d))*newK(d,i);
            fout = newf(i) - newK(i,d)*inv(newK(d,d))*newf(d);
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
            sigmaLockedDofs = [2:3:12, 35:3:obj.nStressDofs]';
            
            sigmaBc = [sigmaLockedDofs, sigmaLockedDofs*0];
            
            sigmaFreeDofs = setdiff(allSigmaDofs,sigmaLockedDofs);
            
            beta = zeros(obj.nStressDofs,1);
            beta(sigmaLockedDofs) = sigmaBc(:,2);
            
            tempF = -(Ae*a + Be*alpha);
            beta(sigmaFreeDofs) = Ce(sigmaFreeDofs, sigmaFreeDofs)\(tempF(sigmaFreeDofs) - Ce(sigmaFreeDofs,sigmaLockedDofs)*sigmaBc(:,2));
        end
        
        % TODO Old method - needs rewriting
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
        
        function dofs = getStressComponentDofs(obj, component, nodes)
            
            ss = obj.stressInterpOrder;
            
            startFrom = 0;
            for i=1:(component-1)
                temp = ss(i)*2*2 ;
                startFrom = startFrom + temp;
            end
            
            dofs = startFrom + nodes;
            
        end
        
        function [bottomDofs, topDofs] = getStressComponentTopAndBottonDofs(obj, component)
            
            ss = obj.stressInterpOrder;
            
            startFrom = 0;
            for i=1:(component-1)
                temp = ss(i)*2*2 ;
                startFrom = startFrom + temp;
            end
            
            temp = 1:obj.stressInterpOrder(component)*2*2;
            
            dofs = startFrom + [temp(1:4), temp((end-3):end) ];
            bottomDofs = dofs(1:4);
            topDofs    = dofs(5:8);
        end
        
        function P = getStressPmatrix(obj, lcoords)
            
            P = [];
            for is = 1:length(obj.stressInterp)
                NvecSigma = obj.stressInterp{is}.eval_N(lcoords);
                P = blkdiag(P,NvecSigma);
            end
            
        end
    end
    
end

