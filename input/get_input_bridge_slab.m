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

% Normal (non-bridge) track sections before/after bridge [m].
in_data.bridge.normal_length_left = 12;
in_data.bridge.normal_length_right = 15;

% Rail fastening spacing [m] for spring layers.
in_data.bridge.fastening_spacing = 0.6;

% Bearing material ID mapping.
% Scalar: same material for all decks.
% Vector: one material ID per deck.
in_data.bridge.bearing_mater_id = 11;

% Slab-ground support material in normal track sections.
in_data.bridge.normal_support_mater_id = 12;

% Vertical offset of bearing ground nodes below bridge nodes [m].
in_data.bridge.bearing_ground_drop = 0.5;

%%
% MATERIAL SPRING DATA: BRIDGE BEARING (for mesh_bridge_slab)
% [K_bearing; C_bearing]
in_data.mater(11).Data = [8.0e8; 1.0e5];
in_data.mater(11).Note = 'bridge bearing';

% MATERIAL SPRING DATA: NORMAL TRACK SLAB-GROUND SUPPORT
% [K_support; C_support]
in_data.mater(12).Data = [5.0e7; 8.0e4];
in_data.mater(12).Note = 'normal section slab-ground support';

end
