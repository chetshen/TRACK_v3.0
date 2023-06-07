%%
%form track matrix in full track model
%%
function [mat] = form_mat_trk_2(inp,geo)
%%
% Initial vectors and matrices with zeros
GDof=2*length(geo.ND);
stiffness=zeros(GDof,GDof);
mass=zeros(GDof,GDof);
damping=zeros(GDof,GDof);
%%
% Caculate the element matrices for each element and assemble them into the
% 
for i = 1:length(geo.EL)
    indice=geo.EL(i,2:3);
    elemType=geo.EL(i,6);
    materID=geo.EL(i,5);
    mater=inp.mater(materID).Data;
    elemLength=node_distance(indice(1),indice(2),geo.ND);
    [Me,Ke,De,elementDof]=form_elem_mat(indice,elemType,elemLength,mater);
    stiffness(elementDof,elementDof)=stiffness(elementDof,elementDof)+Ke;
    mass(elementDof,elementDof)=mass(elementDof,elementDof)+Me;
    damping(elementDof,elementDof)=damping(elementDof,elementDof)+De;
end
%% Boundary conditions
% boundary conditions
fixedNodeU=geo.fixedNodeU;% vertical
fixedNodeV=geo.fixedNodeV; %; rotational

% Active degrees of fredoom
fixedDof=[2*(fixedNodeU-1)+1;2*(fixedNodeV-1)+2];
activeDof=setdiff((1:GDof)',fixedDof);
% col_actDof=length(activeDof);
%%
%STIFFNESS MATRIX REDUCED: considered only the active nodes
mat.K_reduced=sparse(stiffness(activeDof,activeDof));

%MASS MATRIX REDUCED: considered only the active nodes
mat.M_reduced=sparse(mass(activeDof,activeDof));

mat.C_reduced=sparse(damping(activeDof,activeDof));
mat.fixedDof=fixedDof;
mat.activeDof=activeDof;

end

