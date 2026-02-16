%%
% Input file for slab-track bridge model with multiple bridge decks.
% It extends the baseline track input by adding bridge/deck/bearing fields
% consumed by mesh_bridge_slab.m
%%
function [in_data] = get_input_bridge_slab(in_data)

% Start from baseline parameters
in_data = get_input_4(in_data);

%%
% BRIDGE-SLAB MODEL INPUTS (used by mesh_bridge_slab)
% Multiple deck lengths [m]; can have different lengths.
in_data.bridge.deck_lengths = [10, 12, 8];

% Rail fastening spacing [m] for spring layers.
in_data.bridge.fastening_spacing = 0.6;

% Bearing material ID mapping.
% Scalar: same material for all decks.
% Vector: one material ID per deck.
in_data.bridge.bearing_mater_id = 11;

%%
% MATERIAL SPRING DATA: BRIDGE BEARING (for mesh_bridge_slab)
% [K_bearing; C_bearing]
in_data.mater(11).Data = [8.0e8; 1.0e5];
in_data.mater(11).Note = 'bridge bearing';

end
