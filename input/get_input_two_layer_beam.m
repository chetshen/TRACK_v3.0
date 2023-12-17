%%
%Parameter to match measurement in Maider's thesis 20180312
%%
function [in_data,NNslpf] = get_input_two_layer_beam(in_data)
%% Parameters that are relevant to the problem

% Geometry and mesh
in_data.geo.Ltot_R=60;                  % track length[m] 
in_data.geo.SlpSpc=0.6;           % sleeper spacing [m] 
in_data.mesh.numElem_R_betwSprings=6;   %Number of elements between 2 sleepers

% MATERIAL DATA 

in_data.mater(1).Data=[210e9;    %E_R[N/m^2]
                    0.23896e-4;     %I_R[m^4]
                    7.0515e-3;   %A_R[m^2]
                    7800;    %rho_R[kg/m^3]
                      8.1e10;    %G_R
                     0.4];% 0.34] ;     %kappa_R
in_data.mater(1).Note='rail';                  

%MATERIAL SLEEPERS DATA 
in_data.mater(2).Data=[74.6e12; % E  [N/m^2]
                    1.375e-4;  % I [m^4]
                    0.043; % A [m^2]
                    2500;% rho [kg/m^3]
                    74.6e12/2.34; % G
                    0.833]; %0.833
in_data.mater(2).Note='sleeper';

%MATERIAL SPRING-DASHPOT DATA: RAILPAD
in_data.mater(3).Data=[6.5e8; % Stiffness [N/m]
                      6.75e1];% Damping [N.s/m]
in_data.mater(3).Note='railpad';


%MATERIAL SPRING-DASHPOT DATA: BALLAST

in_data.mater(4).Data = [11e7/11;  % Stiffness [N/m]
                         4.5e4/11];% Damping[N.s/m]
in_data.mater(4).Note='ballast';

% solver setting
in_data.solver.n_ts=5000;                      %Number of time steps
in_data.solver.deltat=1/10000;                   % Time step length
in_data.solver.Vx = 30;             % vehicle speed [m/s]
in_data.solver.xw0 = 15;          % initial position of the wheel
%% parameters that are not relevant to the problem
%Geometry input
in_data.geo.dist_S=in_data.geo.SlpSpc; %[m] 
in_data.geo.LExt_S=0.54;             %[m]
in_data.geo.LInt_S=0.750;              %[m]
in_data.geo.TrackWidth=1.5;            %[m]
% in_data.geo.Ltot_S=in_data.geo.LExt_S*2+in_data.geo.TrackWidth;      %[m] full length=2.36[m])
in_data.geo.dist_RS=0.2;               %[m] Between rails and sleepers
in_data.geo.dist_SB=0.2;               %[m] Between sleepers and ballast
in_data.geo.Rw=0.46;                   %[m] Wheel radius
in_data.geo.irr=0;                     % Irregularities
in_data.geo.wb = 2.9;                  %[m] Wheelbase length

%%
%External force

in_data.ext_force.timeh='example.txt';%'FW_h30w40'; %['white_noise.txt']; %[ 'example.txt' ];        %time history of external force
in_data.ext_force.sf=25000;
in_data.ext_force.x=[6,-0.75,0];
in_data.ext_force.Vx=0;
% zdd=load(in_data.ext_force.timeh);
% dof=299;
% in_data.ext_force.dof=ones(length(zdd),1)*dof;       %time history of applied DOF 
in_data.ext_force.wh_ld = 8000*9.8 ; % 8000*9.8 ;% 12742*9.8;   %[N]

%%
%Solver settings

in_data.solver.linsolver_id=2;             %linear solver id, 1 for LDL, 2 for mldivide

%%
%MESH PARAMETERS
% in_data.mesh.numElem_R_betwSprings_L=60;   %Number of elements between 2 springs
in_data.mesh.RefinedMeshLength=0.001;    %Element length at refined mesh around irregularity [m]
in_data.mesh.m_1S_Ext=2;                %Number of elements in a sleeper external
in_data.mesh.m_1S_Int=6;                %Number of elements in a sleeper internal
NEslph=(2*in_data.mesh.m_1S_Ext+in_data.mesh.m_1S_Int)/2; %number of elements for half sleeper
NNslpf=NEslph*2+1;
NNslph=NEslph+1;

in_data.mesh.btypr=1;                   %mesh beam type: 1 for Euler, 2 for Timoshenko
in_data.mesh.btyps=1;
%% Material



%CONTACT
% in_data.mater(5).ElemType=5;   %5 for contact
in_data.mater(5).Data = 8.4e10;%8.4e10;% 8.7e10; %C_Hertz %2/3*in_data.mater.E_R/(1-0.27^2)*sqrt(in_data.geo.Rw) ;% 8.7e10; %[N/m^(3/2)] from (Steffens, 2005)
in_data.mater(5).Note='contact';
%VEHICLE
% in_data.mater(6).ElemType=6;   %6 for vehicle
in_data.mater(6).Data =  [0;   %M_sprg%4653.5[kg]
                      883.6;%883.6; %M_unsprg[kg]
                     1.22e6]; %K_PS[N/m]
in_data.mater(6).wsfile=[];%'E:\FE\Model 1 wheel\half_ws.spm';  %.spm file generated from ANSYS; empty means a rigid wheelset               
in_data.mater(6).Note='vehicle';
                 
%==================Add new materials below=================================
% in_data.mater(7).ElemType=2;       %1 for rail material; 2 for sleeper; 3 for railpads; 4 for ballast
in_data.mater(7).Data=[210e9;    %E_R[N/m^2]
                   2.3379e-5;    %I_R[m^4]
                    6.977e-3;    %A_R[m^3]
                        7800;    %rho_R[kg/m^3]
                      8.1e10;    %G_R
                      0.34] ;     %kappa_R
in_data.mater(7).Note='rail degraded';
%CONTACT
in_data.mater(8).Data = 2/3*in_data.mater(1).Data(1)/(1-0.27^2)*sqrt(0.3) ;% 8.7e10; %[N/m^(3/2)] from (Steffens, 2005)
in_data.mater(8).Note='contact parameter used for winkler bedding';

in_data.mater(9).Data=[1.100000000000000e+9;0;0;8.836000000000000e+02];
in_data.mater(9).Note='mass spring model';

% 
in_data.mater(10).Data = 1.1e9;%1.1e9; %C_Hertz %2/3*in_data.mater.E_R/(1-0.27^2)*sqrt(in_data.geo.Rw) ;% 8.7e10; %[N/m^(3/2)] from (Steffens, 2005)
in_data.mater(10).Note='linear contact';


end