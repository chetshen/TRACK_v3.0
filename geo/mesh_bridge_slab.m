function [geo,nodeCoord] = mesh_bridge_slab(in_data,beam_type_rail,beam_type_slab,beam_type_bridge)
%%
% Build a slab-track bridge model with layered components (top to bottom):
%   1) rail beam
%   2) rail fastening spring/damper layer (0.6 m spacing by default)
%   3) slab beam
%   4) slab-bridge interlayer spring/damper layer
%   5) bridge deck beam (multiple decks/spans supported by bearings)
%
% Multiple bridge decks can be defined in input as:
%   in_data.bridge.deck_lengths = [L1 L2 ... Ln];  % [m]
% Each deck is supported at both ends by bearing springs/dampers.
% Bearing material(s) are referenced from in_data.mater via:
%   in_data.bridge.bearing_mater_id = scalar or [id1 ... idn]
%
% Optional settings:
%   in_data.bridge.fastening_spacing           (default 0.6 m)
%   in_data.mesh.numElem_R_betwSprings         (default 1)
%
% Material IDs used in geo.EL(:,5):
%   1 rail beam
%   2 slab beam
%   3 bridge deck beam
%   4 rail fastenings (spring/damper)
%   5 slab-bridge interlayer (spring/damper)
%   6 bridge bearings (spring/damper, material ID from in_data.bridge.bearing_mater_id)
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

% Bearing-material mapping for bearings (defined in get_input_*.m)
% Use in_data.bridge.bearing_mater_id as scalar or one value per deck.
if ~isfield(in_data,'bridge') || ~isfield(in_data.bridge,'bearing_mater_id')
    error('Define in_data.bridge.bearing_mater_id in get_input_*.m.');
end
bearing_mater_id = expand_to_num_decks(in_data.bridge.bearing_mater_id, numel(deck_lengths), 'bearing_mater_id');

%% Build x coordinates
Ltot = sum(deck_lengths);
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

x_support = [0, cumsum(deck_lengths)];
x_support = unique(round(x_support,10));

% Ensure support locations always exist in beam mesh
x_beam = unique(sort([x_beam, x_support]));

nBeam = numel(x_beam);

%% Nodes: rail/slab/bridge layers + support-ground nodes for bearings
nodeCoord_R = [x_beam(:), zeros(nBeam,1), zeros(nBeam,1), ones(nBeam,1)];
nodeCoord_S = [x_beam(:), zeros(nBeam,1), -dist_RS*ones(nBeam,1), 2*ones(nBeam,1)];
nodeCoord_B = [x_beam(:), zeros(nBeam,1), -(dist_RS+dist_SB)*ones(nBeam,1), 3*ones(nBeam,1)];

% Ground/support nodes for bearing springs
nSupport = numel(x_support);
nodeCoord_G = [x_support(:), zeros(nSupport,1), -(dist_RS+dist_SB)*ones(nSupport,1), 4*ones(nSupport,1)];

nodeCoord = [nodeCoord_R;nodeCoord_S;nodeCoord_B;nodeCoord_G];
geo.ND = [(1:size(nodeCoord,1))', nodeCoord];

%% Elements
% [node1 node2 partID materialID elemType]
elemNodes = zeros(0,5);

% Beam layers
railBeam = [(1:nBeam-1)', (2:nBeam)', ones(nBeam-1,1), ones(nBeam-1,1), beam_type_rail*ones(nBeam-1,1)];
slabStart = nBeam;
slabBeam = [slabStart+(1:nBeam-1)', slabStart+(2:nBeam)', 2*ones(nBeam-1,1), 2*ones(nBeam-1,1), beam_type_slab*ones(nBeam-1,1)];
bridgeStart = 2*nBeam;
bridgeBeam = [bridgeStart+(1:nBeam-1)', bridgeStart+(2:nBeam)', 3*ones(nBeam-1,1), 3*ones(nBeam-1,1), beam_type_bridge*ones(nBeam-1,1)];

elemNodes = [elemNodes; railBeam; slabBeam; bridgeBeam]; %#ok<AGROW>

% Rail fastenings (rail-slab)
[xFoundSpr, railNodeIdx] = ismember(round(x_layer_spr,10), round(x_beam,10));
if ~all(xFoundSpr)
    error('Unable to map fastening locations to beam nodes.');
end
slabNodeIdx = nBeam + railNodeIdx;
fasteningElem = [railNodeIdx(:), slabNodeIdx(:), 4*ones(numel(railNodeIdx),1), 4*ones(numel(railNodeIdx),1), 3*ones(numel(railNodeIdx),1)];
elemNodes = [elemNodes; fasteningElem]; %#ok<AGROW>

% Slab-bridge interlayer
bridgeNodeIdx = 2*nBeam + railNodeIdx;
interlayerElem = [slabNodeIdx(:), bridgeNodeIdx(:), 5*ones(numel(railNodeIdx),1), 5*ones(numel(railNodeIdx),1), 3*ones(numel(railNodeIdx),1)];
elemNodes = [elemNodes; interlayerElem]; %#ok<AGROW>

% Bearings at deck ends (two bearings per deck)
deckEndA = [0, cumsum(deck_lengths(1:end-1))];
deckEndB = cumsum(deck_lengths);

[xFoundA, idxA] = ismember(round(deckEndA,10), round(x_beam,10));
[xFoundB, idxB] = ismember(round(deckEndB,10), round(x_beam,10));
idxG = 1:nSupport;
if ~all(xFoundA) || ~all(xFoundB)
    error('Unable to map bearing locations to node coordinates.');
end

% Ground/support node ids in global indexing
groundStart = 3*nBeam;
groundNodeIdx = groundStart + idxG;

bearingElem = zeros(2*numel(deck_lengths),5);
for iDeck = 1:numel(deck_lengths)
    deckNodeA = 2*nBeam + idxA(iDeck);
    deckNodeB = 2*nBeam + idxB(iDeck);

    gNodeA = groundNodeIdx(find(abs(x_support - deckEndA(iDeck)) < 1e-12,1,'first'));
    gNodeB = groundNodeIdx(find(abs(x_support - deckEndB(iDeck)) < 1e-12,1,'first'));

    bearingElem(2*iDeck-1,:) = [deckNodeA, gNodeA, 6, bearing_mater_id(iDeck), 3];
    bearingElem(2*iDeck,:)   = [deckNodeB, gNodeB, 6, bearing_mater_id(iDeck), 3];
end

elemNodes = [elemNodes; bearingElem]; %#ok<AGROW>

geo.EL = [(1:size(elemNodes,1))', elemNodes];

%% Element counters
m_RailBeam = size(railBeam,1);
m_SlabBeam = size(slabBeam,1);
m_BridgeBeam = size(bridgeBeam,1);
m_Fastening = size(fasteningElem,1);
m_Interlayer = size(interlayerElem,1);
m_Bearing = size(bearingElem,1);
geo.NumEL = [m_RailBeam,m_SlabBeam,m_BridgeBeam,m_Fastening,m_Interlayer,m_Bearing,size(elemNodes,1)];

%% Boundary conditions
% - Support-ground nodes are fixed in U and V (anchor for bearing springs)
% - Rail end nodes are fixed in U and V to prevent rigid-body drift
railNodes = 1:nBeam;
supportNodes = (3*nBeam+1):(3*nBeam+nSupport);

geo.fixedNodeU = [railNodes(1); railNodes(end); supportNodes(:)];
geo.fixedNodeV = [railNodes(1); railNodes(end); supportNodes(:)];

%% Useful indexing
geo.layerNodes.rail = (1:nBeam)';
geo.layerNodes.slab = (nBeam+1:2*nBeam)';
geo.layerNodes.bridge = (2*nBeam+1:3*nBeam)';
geo.layerNodes.support = (3*nBeam+1:3*nBeam+nSupport)';
geo.deck.supportX = x_support(:);
geo.deck.lengths = deck_lengths(:);
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
