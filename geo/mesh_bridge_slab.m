function [geo,nodeCoord] = mesh_bridge_slab(in_data,beam_type_rail,beam_type_slab,beam_type_bridge)
%%
% Build a slab-track bridge model with layered components (top to bottom):
%   1) rail beam
%   2) rail fastening spring/damper layer (0.6 m spacing by default)
%   3) slab beam
%   4) bridge zone: slab-bridge interlayer spring/damper layer
%      approach zones: slab-ground spring/damper layer
%   5) bridge deck beam (multiple decks/spans supported by bearings)
%
% Multiple bridge decks can be defined in input as:
%   in_data.bridge.deck_lengths = [L1 L2 ... Ln];  % [m]
% Approach track sections:
%   in_data.bridge.normal_length_left  = L_left;   % [m]
%   in_data.bridge.normal_length_right = L_right;  % [m]
%
% Optional settings:
%   in_data.bridge.fastening_spacing           (default 0.6 m)
%   in_data.mesh.numElem_R_betwSprings         (default 1)
%   in_data.bridge.normal_support_mater_id     (default 5)
%
% Material IDs used in geo.EL(:,5):
%   1 rail beam
%   2 slab beam
%   3 bridge deck beam
%   4 rail fastenings (spring/damper)
%   5 slab-bridge interlayer (spring/damper)
%   6 bridge bearings (spring/damper, material ID from in_data.bridge.bearing_mater_id)
%   7 slab-ground support in normal sections (material ID from in_data.bridge.normal_support_mater_id)
%%

if nargin < 2 || isempty(beam_type_rail)
    beam_type_rail = 2;
end
if nargin < 3 || isempty(beam_type_slab)
    beam_type_slab = beam_type_rail;
end
if nargin < 4 || isempty(beam_type_bridge)
    beam_type_bridge = beam_type_slab;
end

% Geometry / mesh inputs
if ~isfield(in_data,'bridge') || ~isfield(in_data.bridge,'deck_lengths')
    if isfield(in_data,'geo') && isfield(in_data.geo,'Ltot_R')
        deck_lengths = in_data.geo.Ltot_R;
    else
        error('Define in_data.bridge.deck_lengths (or in_data.geo.Ltot_R for a single deck).');
    end
else
    deck_lengths = in_data.bridge.deck_lengths;
end

deck_lengths = deck_lengths(:)';
if isempty(deck_lengths) || any(deck_lengths <= 0)
    error('in_data.bridge.deck_lengths must contain positive values.');
end

if isfield(in_data,'bridge') && isfield(in_data.bridge,'normal_length_left')
    normal_length_left = in_data.bridge.normal_length_left;
else
    normal_length_left = 0;
end
if isfield(in_data,'bridge') && isfield(in_data.bridge,'normal_length_right')
    normal_length_right = in_data.bridge.normal_length_right;
else
    normal_length_right = 0;
end
if normal_length_left < 0 || normal_length_right < 0
    error('in_data.bridge.normal_length_left/right must be nonnegative.');
end

if isfield(in_data,'bridge') && isfield(in_data.bridge,'fastening_spacing')
    fastening_spacing = in_data.bridge.fastening_spacing;
elseif isfield(in_data,'geo') && isfield(in_data.geo,'SlpSpc')
    fastening_spacing = in_data.geo.SlpSpc;
else
    fastening_spacing = 0.6;
end

if fastening_spacing <= 0
    error('Fastening spacing must be positive.');
end

if isfield(in_data,'mesh') && isfield(in_data.mesh,'numElem_R_betwSprings')
    numElem_between_fastenings = in_data.mesh.numElem_R_betwSprings;
else
    numElem_between_fastenings = 1;
end
if numElem_between_fastenings < 1 || mod(numElem_between_fastenings,1) ~= 0
    error('numElem_R_betwSprings must be a positive integer.');
end

% Vertical distances
if ~isfield(in_data,'geo') || ~isfield(in_data.geo,'dist_RS') || ~isfield(in_data.geo,'dist_SB')
    error('in_data.geo.dist_RS and in_data.geo.dist_SB are required.');
end
dist_RS = in_data.geo.dist_RS;
dist_SB = in_data.geo.dist_SB;

% Bearing material mapping
if ~isfield(in_data,'bridge') || ~isfield(in_data.bridge,'bearing_mater_id')
    error('Define in_data.bridge.bearing_mater_id in get_input_*.m.');
end
bearing_mater_id = expand_to_num_decks(in_data.bridge.bearing_mater_id, numel(deck_lengths), 'bearing_mater_id');

% Normal slab support material mapping
if isfield(in_data,'bridge') && isfield(in_data.bridge,'normal_support_mater_id')
    normal_support_mater_id = in_data.bridge.normal_support_mater_id;
else
    normal_support_mater_id = 5;
end

%% Build x coordinates
Lbridge = sum(deck_lengths);
xBridgeStart = normal_length_left;
xBridgeEnd = normal_length_left + Lbridge;
Ltot = normal_length_left + Lbridge + normal_length_right;

dx = fastening_spacing / numElem_between_fastenings;

x_beam = 0:dx:Ltot;
if abs(x_beam(end)-Ltot) > 1e-12
    x_beam = [x_beam, Ltot];
end
x_beam = unique(round(x_beam,10));

x_layer_spr = 0:fastening_spacing:Ltot;
if abs(x_layer_spr(end)-Ltot) > 1e-12
    x_layer_spr = [x_layer_spr, Ltot];
end
x_layer_spr = unique(round(x_layer_spr,10));

x_support = xBridgeStart + [0, cumsum(deck_lengths)];
x_support = unique(round(x_support,10));

% Ensure special positions exist in the rail/slab mesh
x_beam = unique(sort([x_beam, x_support, xBridgeStart, xBridgeEnd]));

nBeam = numel(x_beam);
nDeck = numel(deck_lengths);

%% Nodes: rail/slab/bridge + bearing ground + normal-support ground
nodeCoord_R = [x_beam(:), zeros(nBeam,1), zeros(nBeam,1), ones(nBeam,1)];
nodeCoord_S = [x_beam(:), zeros(nBeam,1), -dist_RS*ones(nBeam,1), 2*ones(nBeam,1)];

% Bridge nodes are built deck-by-deck so adjacent decks do NOT share nodes.
deckStartX = xBridgeStart + [0, cumsum(deck_lengths(1:end-1))];
deckEndX = xBridgeStart + cumsum(deck_lengths);

nodeCoord_B = zeros(0,4);
bridgeElem = zeros(0,5);
bridgeNodeDeckStart = zeros(1,nDeck);
bridgeNodeDeckEnd = zeros(1,nDeck);

baseBridgeNode = 2*nBeam;
for iDeck = 1:nDeck
    xA = round(deckStartX(iDeck),10);
    xB = round(deckEndX(iDeck),10);

    idxA = find(abs(x_beam - xA) < 1e-12,1,'first');
    idxB = find(abs(x_beam - xB) < 1e-12,1,'first');
    if isempty(idxA) || isempty(idxB)
        error('Unable to map deck boundaries to beam coordinates.');
    end

    xDeck = x_beam(idxA:idxB);
    nDeckNode = numel(xDeck);

    thisStart = baseBridgeNode + size(nodeCoord_B,1) + 1;
    thisNodeIds = thisStart:(thisStart+nDeckNode-1);

    nodeCoord_B = [nodeCoord_B; ...
        xDeck(:), zeros(nDeckNode,1), -(dist_RS+dist_SB)*ones(nDeckNode,1), 3*ones(nDeckNode,1)]; %#ok<AGROW>

    if nDeckNode >= 2
        bridgeElem = [bridgeElem; ...
            thisNodeIds(1:end-1)', thisNodeIds(2:end)', 3*ones(nDeckNode-1,1), 3*ones(nDeckNode-1,1), beam_type_bridge*ones(nDeckNode-1,1)]; %#ok<AGROW>
    end

    bridgeNodeDeckStart(iDeck) = thisNodeIds(1);
    bridgeNodeDeckEnd(iDeck) = thisNodeIds(end);
end
nBridge = size(nodeCoord_B,1);

% Ground nodes for bearing springs
nSupport = numel(x_support);
nodeCoord_Gb = [x_support(:), zeros(nSupport,1), -(dist_RS+dist_SB)*ones(nSupport,1), 4*ones(nSupport,1)];

% Ground nodes for normal-section slab supports
isNormalSpr = (x_layer_spr < xBridgeStart-1e-12) | (x_layer_spr > xBridgeEnd+1e-12);
x_normal_spr = x_layer_spr(isNormalSpr);
nNormalSpr = numel(x_normal_spr);
nodeCoord_Gn = [x_normal_spr(:), zeros(nNormalSpr,1), -dist_RS*ones(nNormalSpr,1), 5*ones(nNormalSpr,1)];

nodeCoord = [nodeCoord_R;nodeCoord_S;nodeCoord_B;nodeCoord_Gb;nodeCoord_Gn];
geo.ND = [(1:size(nodeCoord,1))', nodeCoord];

%% Elements
% [node1 node2 partID materialID elemType]
elemNodes = zeros(0,5);

% Rail/slab beam layers
railBeam = [(1:nBeam-1)', (2:nBeam)', ones(nBeam-1,1), ones(nBeam-1,1), beam_type_rail*ones(nBeam-1,1)];
slabStart = nBeam;
slabBeam = [slabStart+(1:nBeam-1)', slabStart+(2:nBeam)', 2*ones(nBeam-1,1), 2*ones(nBeam-1,1), beam_type_slab*ones(nBeam-1,1)];

elemNodes = [elemNodes; railBeam; slabBeam; bridgeElem]; %#ok<AGROW>

% Rail fastenings (rail-slab)
[xFoundSpr, railNodeIdx] = ismember(round(x_layer_spr,10), round(x_beam,10));
if ~all(xFoundSpr)
    error('Unable to map fastening locations to beam nodes.');
end
slabNodeIdx = nBeam + railNodeIdx;
fasteningElem = [railNodeIdx(:), slabNodeIdx(:), 4*ones(numel(railNodeIdx),1), 4*ones(numel(railNodeIdx),1), 3*ones(numel(railNodeIdx),1)];
elemNodes = [elemNodes; fasteningElem]; %#ok<AGROW>

% Bridge-zone slab-bridge interlayer springs
isBridgeSpr = ~isNormalSpr;
x_bridge_spr = round(x_layer_spr(isBridgeSpr),10);
interlayerElem = zeros(0,5);
for iDeck = 1:nDeck
    xA = round(deckStartX(iDeck),10);
    xB = round(deckEndX(iDeck),10);

    inDeck = (x_bridge_spr >= xA-1e-12) & (x_bridge_spr <= xB+1e-12);
    if iDeck > 1
        inDeck = inDeck & (x_bridge_spr > xA+1e-12); % avoid duplicating shared deck boundary
    end

    xDeckSpr = x_bridge_spr(inDeck);
    [okSlab, slabIdxDeck] = ismember(xDeckSpr, round(x_beam,10));
    if ~all(okSlab)
        error('Unable to map interlayer spring positions to slab nodes.');
    end

    bridgeNodeRange = bridgeNodeDeckStart(iDeck):bridgeNodeDeckEnd(iDeck);
    xBridgeDeck = round(nodeCoord_B(bridgeNodeRange-baseBridgeNode,1)',10);
    [okBridge, locBridge] = ismember(xDeckSpr, xBridgeDeck);
    if ~all(okBridge)
        error('Unable to map interlayer spring positions to bridge nodes.');
    end

    bridgeIdxDeck = bridgeNodeRange(locBridge);
    slabIdxDeck = nBeam + slabIdxDeck;

    interlayerElem = [interlayerElem; ...
        slabIdxDeck(:), bridgeIdxDeck(:), 5*ones(numel(xDeckSpr),1), 5*ones(numel(xDeckSpr),1), 3*ones(numel(xDeckSpr),1)]; %#ok<AGROW>
end
elemNodes = [elemNodes; interlayerElem]; %#ok<AGROW>

% Bearings at deck ends (two bearings per deck)
groundStartBearing = 2*nBeam + nBridge;
groundNodeBearing = groundStartBearing + (1:nSupport);

bearingElem = zeros(2*nDeck,5);
for iDeck = 1:nDeck
    gNodeA = groundNodeBearing(find(abs(x_support - deckStartX(iDeck)) < 1e-12,1,'first'));
    gNodeB = groundNodeBearing(find(abs(x_support - deckEndX(iDeck)) < 1e-12,1,'first'));

    bearingElem(2*iDeck-1,:) = [bridgeNodeDeckStart(iDeck), gNodeA, 6, bearing_mater_id(iDeck), 3];
    bearingElem(2*iDeck,:)   = [bridgeNodeDeckEnd(iDeck),   gNodeB, 6, bearing_mater_id(iDeck), 3];
end
elemNodes = [elemNodes; bearingElem]; %#ok<AGROW>

% Normal-section slab-ground support springs
groundStartNormal = groundStartBearing + nSupport;
groundNodeNormal = groundStartNormal + (1:nNormalSpr);
normalSupportElem = zeros(nNormalSpr,5);
if nNormalSpr > 0
    [okNormalSlab, slabIdxNormal] = ismember(round(x_normal_spr,10), round(x_beam,10));
    if ~all(okNormalSlab)
        error('Unable to map normal-section spring positions to slab nodes.');
    end
    slabIdxNormal = nBeam + slabIdxNormal;

    normalSupportElem = [slabIdxNormal(:), groundNodeNormal(:), 7*ones(nNormalSpr,1), normal_support_mater_id*ones(nNormalSpr,1), 3*ones(nNormalSpr,1)];
    elemNodes = [elemNodes; normalSupportElem]; %#ok<AGROW>
end

geo.EL = [(1:size(elemNodes,1))', elemNodes];

%% Element counters
m_RailBeam = size(railBeam,1);
m_SlabBeam = size(slabBeam,1);
m_BridgeBeam = size(bridgeElem,1);
m_Fastening = size(fasteningElem,1);
m_Interlayer = size(interlayerElem,1);
m_Bearing = size(bearingElem,1);
m_NormalSupport = size(normalSupportElem,1);
geo.NumEL = [m_RailBeam,m_SlabBeam,m_BridgeBeam,m_Fastening,m_Interlayer,m_Bearing,m_NormalSupport,size(elemNodes,1)];

%% Boundary conditions
% - Bearing-ground nodes and normal-ground nodes are fixed in U and V
% - Rail end nodes are fixed in U and V to prevent rigid-body drift
railNodes = 1:nBeam;
supportNodesBearing = (groundStartBearing+1):(groundStartBearing+nSupport);
supportNodesNormal = (groundStartNormal+1):(groundStartNormal+nNormalSpr);

geo.fixedNodeU = [railNodes(1); railNodes(end); supportNodesBearing(:); supportNodesNormal(:)];
geo.fixedNodeV = [railNodes(1); railNodes(end); supportNodesBearing(:); supportNodesNormal(:)];

%% Useful indexing
geo.layerNodes.rail = (1:nBeam)';
geo.layerNodes.slab = (nBeam+1:2*nBeam)';
geo.layerNodes.bridge = (2*nBeam+1:2*nBeam+nBridge)';
geo.layerNodes.supportBearing = (groundStartBearing+1:groundStartBearing+nSupport)';
geo.layerNodes.supportNormal = (groundStartNormal+1:groundStartNormal+nNormalSpr)';
geo.deck.supportX = x_support(:);
geo.deck.lengths = deck_lengths(:);
geo.deck.bridgeNodeStart = bridgeNodeDeckStart(:);
geo.deck.bridgeNodeEnd = bridgeNodeDeckEnd(:);
geo.bridgeStartX = xBridgeStart;
geo.bridgeEndX = xBridgeEnd;
geo.normalLengthLeft = normal_length_left;
geo.normalLengthRight = normal_length_right;
geo.fasteningSpacing = fastening_spacing;

end

function vec = expand_to_num_decks(value,numDecks,name)
vec = value(:)';
if isscalar(vec)
    vec = repmat(vec,1,numDecks);
elseif numel(vec) ~= numDecks
    error('in_data.bridge.%s must be scalar or have one value per deck.',name);
end
end
