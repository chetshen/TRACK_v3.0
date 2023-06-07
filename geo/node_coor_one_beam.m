function [nodeCoord]=node_coor_one_beam(in_data)
%%
% function to calculate node coordinate based on input data
% Node coordinates for left rail, right rail, sleepers and ballast are
% output in 4 matrix.
% Each matrix includes 4 columns, representing x,y,z-coordinate and
% identifier, respectively.
% The identifiers are defined as
%                1        nodes on left rail
%                2        nodes on right rail
%                3        nodes on sleepers
%                4        nodes on ballast/ground
%%
%Read in data
Ltot = in_data.geo.Ltot_R;                  %[m]
% dist_R_betwSprings=in_data.geo.SlpSpc;      %[m]
numElem=in_data.mesh.numElem_R_betwSprings;   %Number of elements between 2 springs
LInt_S=in_data.geo.TrackWidth;              %[m]

%%
%
% m_Spring_RS=(Ltot/dist_R_betwSprings)+1; %number of springs between ONE rail and sleeper
%

%

%element length
%
% LElem_R=Ltot_R/((numElem_R_betwSprings)*(m_Spring_RS-1)); %element length

LElem_R=Ltot/numElem; %element length

%
%number of elements 
%
m_R=Ltot/LElem_R; %number of elements in ONE RAIL
m_R = int64(m_R);
%Node coordinates for rail
%
nodeCoord = [(0:LElem_R:Ltot);-LInt_S/2.*ones(1,(m_R+1));zeros(1,(m_R+1));ones(1,(m_R+1))]';

