%%
% Mesh one beam
function [geo] = mesh_one_beam(beam_type_r,nodeCoord)

[nodeCoord_beam,~,~,~]=retrieve_coord(nodeCoord);
%[nodeCoord_Rl,nodeCoord_Rr,nodeCoord_S,nodeCoord_B,nodeCoord]=node_coor(in_data);
geo.ND = [1:length(nodeCoord);nodeCoord']';

%%
%number of elements for each part
m_beam=length(nodeCoord_beam)-1;
m_total = m_beam;

%%
% elements and corresponding
% [element number, node1,node2, partID, materialID, element type]
%Element type:
%          1     Euler beam
%          2     Timo beam
%          3     Spring/damper pair
%          4     Contact
%          5     Vehicle

elemNodes = zeros(m_total,5);

for i=1:m_beam %left Rail
    elemNodes(i,:)=[i,i+1,1,1,beam_type_r];
end

%Write to output
%
geo.EL = [1:length(elemNodes);elemNodes']';
geo.NumEL = m_total;

%%
%BOUNDARY CONDITION
%
LeftRailNodes = geo.ND(geo.ND(:,5)==1,1);

geo.fixedNodeU=[LeftRailNodes(1);LeftRailNodes(length(LeftRailNodes))];
%%%%%%%%%%%%%%%
geo.fixedNodeV= [];%[LeftRailNodes(1);LeftRailNodes(length(LeftRailNodes))];
%%%%%%%%%%%%%%%%
end
%%

