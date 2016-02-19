function [f, bc] = cubeBC( opt, f, P, dof, side1nodes, side2nodes, side3nodes, side4nodes, side5nodes)

%P is force
%F is forcevector
%SideXnodes is according to CubeMesher

if strcmp('Konsol', opt)
    side2zdof = dof(side2nodes,3);
    f(side2zdof) = f(side2zdof)+ P/length(side2zdof);

    ldofs = dof(side1nodes,:);
    ldofs = ldofs(:);
    bc = [ldofs,ldofs*0];
elseif strcmp('KonsolMedUtbredd', opt)
    side5zdof = dof(side5nodes,3);
    f(side5zdof) = f(side5zdof) + P/length(side5zdof);

    ldofs = dof(side1nodes,:);
    ldofs = ldofs(:);
    bc = [ldofs,ldofs*0];
elseif(strcmp('FrittUpplagt', opt))
    side5zdof = dof(side5nodes,3);
    f(side5zdof) = P/length(side5zdof);
    
    ldofs2 = dof(side2nodes,1:3);
    ldofs2 = ldofs2(:);
    ldofs = dof(side1nodes,:);
    ldofs = ldofs(:);
    bc = [[ldofs;ldofs2],[ldofs;ldofs2]*0];
elseif(strcmp('InspandPlatta', opt))
    side5zdof = dof(side5nodes,3);
    f(side5zdof) = f(side5zdof) +  P/length(side5zdof);
    
    ldofs = dof([side1nodes;side2nodes;side3nodes;side4nodes], 1:3);
    ldofs = ldofs(:);
    bc = [ldofs, ldofs*0];
end

end

