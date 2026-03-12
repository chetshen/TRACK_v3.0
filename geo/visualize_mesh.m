function h = visualize_mesh(geo, showNodeId, showElemId)
%% visualize_mesh
% Visualize mesh nodes and elements from geo.ND / geo.EL.
%
% Usage:
%   visualize_mesh(geo)
%   visualize_mesh(geo,1,0)
%
% Inputs:
%   geo       : model geometry struct with fields ND and EL
%   showNodeId: 1 to plot node labels (default 0)
%   showElemId: 1 to plot element labels (default 0)

if nargin < 2 || isempty(showNodeId)
    showNodeId = 0;
end
if nargin < 3 || isempty(showElemId)
    showElemId = 0;
end

if ~isfield(geo,'ND') || ~isfield(geo,'EL')
    error('Input geo must contain fields ND and EL.');
end

ND = geo.ND;
EL = geo.EL;

% ND format: [nodeID x y z layer]
% EL format: [elemID node1 node2 partID materialID elemType]

figure;
hold on;
grid on;
box on;
axis equal;

% Plot nodes
scatter3(ND(:,2), ND(:,3), ND(:,4), 16, 'k', 'filled');

if showNodeId == 1
    for i = 1:size(ND,1)
        text(ND(i,2), ND(i,3), ND(i,4), sprintf('N%d',ND(i,1)), ...
            'FontSize', 7, 'Color', [0.1 0.1 0.1]);
    end
end

% Color by element part ID when available
partIds = unique(EL(:,4));
cc = lines(max(numel(partIds),7));

for i = 1:size(EL,1)
    n1 = EL(i,2);
    n2 = EL(i,3);

    p = EL(i,4);
    cIdx = find(partIds == p,1,'first');
    if isempty(cIdx)
        cIdx = 1;
    end

    x = [ND(n1,2), ND(n2,2)];
    y = [ND(n1,3), ND(n2,3)];
    z = [ND(n1,4), ND(n2,4)];
    plot3(x,y,z,'-','Color',cc(cIdx,:),'LineWidth',1.2);

    if showElemId == 1
        xc = 0.5*(x(1)+x(2));
        yc = 0.5*(y(1)+y(2));
        zc = 0.5*(z(1)+z(2));
        text(xc,yc,zc,sprintf('E%d',EL(i,1)), 'FontSize',7, 'Color',[0 0 1]);
    end
end

% Plot fixed nodes if available
if isfield(geo,'fixedNodeU')
    idx = unique(geo.fixedNodeU(:));
    scatter3(ND(idx,2), ND(idx,3), ND(idx,4), 40, 'm', '^');
end
if isfield(geo,'fixedNodeV')
    idx = unique(geo.fixedNodeV(:));
    scatter3(ND(idx,2), ND(idx,3), ND(idx,4), 40, 'g', '>');
end

xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
title('Mesh visualization (nodes + elements)');
view(45,30);

h = gca;
end
