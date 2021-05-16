function [H1, fxx] = black_box_model(kp,cp,kb,cb)
%%
% black_box_model: calculates FRF based on hammer test
% INPUT:
%   kp - railpad stiffness
%   cp - railpad damping
%   kb - ballast stiffness
%   cb - ballast damping
% OUTPUT:
% 	H1 - FRF (point receptance, complex value)
%  fxx - frequency points correponding to FRF
% USAGE EXAMPLE:
% see example_black_box_optimization.mlx
% 

% get input parameters
[inp,NNslpf] = get_input_4(); % default input parameters

% geometry
[nodeCoord] = node_coor(inp); % form nodes
btypr = inp.mesh.btypr;       % element type for rail
btyps = inp.mesh.btyps;       % element type for sleeper
[geo] = mesh_trk_full(btypr,btyps,nodeCoord); % mesh geometry

% modify railpad and ballast stiffness and damping (parameters to be optimized)
inp.mater(3).Data(1) = kp;
inp.mater(3).Data(2) = cp;
inp.mater(4).Data(1) = kb / NNslpf;
inp.mater(4).Data(2) = cb / NNslpf;

% form system matrix
% mat_ws=form_mat_ws(inp);                      % wheelset matrix
mat_trk = form_mat_trk_2(inp,geo);            % form track matrix

% calculate hammer response
[H1,fxx] = hammer_main_snst_mcs(inp,mat_trk,geo);