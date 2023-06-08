%%
%Parameters of the case study on P176 of VehicleBridge Interaction Dynamics with
%applications to high-speed railways by Yang Y.B.
%%
function [in_data,NNslpf] = get_input_moving_load(in_data)
%% Parameters that are relevant to the problem

% Geometry data
in_data.geo.Ltot_R = 25;                    % Bridge length[m]
in_data.mesh.numElem_R_betwSprings = 500;   % Number of elements

% Solver settings
in_data.solver.n_ts = 1800;                       % Number of time steps
in_data.solver.deltat = 1/1000;                   % Time step length [s]
in_data.solver.Vx = (100/3.6);                     % Vehicle speed [m/s]
in_data.solver.xw0 = 0.1;                         % Initial position [m]

% Material data
in_data.mater(1).Data=[2.87e9;                 %E [N/m^2]
                    2.9;                       %I [m^4]
                    1;                         % A [m^2]
                    2303;                      % rho [kg/m^3]
                      2.87e9/(2*(1 + 0.2));    % G
                     0.4];                     % kappa
in_data.mater(1).Note='bridge beam'; 

in_data.mater(7).Data=[2.87e9*0.6;                 %E [N/m^2]
                    2.9;                       %I [m^4]
                    1;                         % A [m^2]
                    2303;                      % rho [kg/m^3]
                      2.87e9/(2*(1 + 0.2));    % G
                     0.4];                     % kappa
in_data.mater(7).Note='damaged beam';

in_data.mater(6).Data =  [0;   % not relevant
                      5750;  % Sprung mass [kg]
                     0]; % not relevant
in_data.mater(6).Note='vehicle';

in_data.mater(10).Data = 1595e3;% Linear spring K
in_data.mater(10).Note='linear contact';
%% Parameters that are not relevant to the problem
%Geometry input


in_data.geo.SlpSpc=0.5;           %[m] 
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
in_data.ext_force.wh_ld = 0; % 8000*9.8 ;% 12742*9.8;   %[N]

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
%%

%MATERIAL  DATA                  

in_data.mater(2).Data=[74.6e9/1.5; % #19.4e12# or #19.4e9# [N/m^2]
                    1.375e-4;  %[m^4]
                    0.043; %[m^2]
                    2500;%2140;%3070;   %2480[kg/m^3]
                    74.6e9/2.34; %E/2.34
                    0.833]; %0.833
in_data.mater(2).Note='sleeper';

in_data.mater(3).Data=[6.5e8; %1.56e9;%1.3e9; %K_Spring_RS [N/m]
                      6.75e4];%6.75e4]; %4.5e4]; %C_Damper_RS[N.s/m]
in_data.mater(3).Note='railpad';


in_data.mater(4).Data = [11e7/NNslpf;  %K_Spring_SB[N/m]
                         4.5e4/NNslpf];%3.444e4/5];  %C_Damper_SB[N.s/m]
in_data.mater(4).Note='ballast';


%CONTACT
in_data.mater(5).Data = 8.4e10; %K_Hertz %2/3*in_data.mater.E_R/(1-0.27^2)*sqrt(in_data.geo.Rw) ;% 8.7e10; %[N/m^(3/2)] from (Steffens, 2005)
in_data.mater(5).Note='contact';



in_data.mater(6).wsfile=[];%'E:\FE\Model 1 wheel\half_ws.spm';  %.spm file generated from ANSYS; empty means a rigid wheelset               

                 
%==================Add new materials below=================================

%CONTACT
in_data.mater(8).Data = 2/3*in_data.mater(1).Data(1)/(1-0.27^2)*sqrt(0.3) ;% 8.7e10; %[N/m^(3/2)] from (Steffens, 2005)
in_data.mater(8).Note='contact parameter used for winkler bedding';

in_data.mater(9).Data=[1.100000000000000e+9;0;0;8.836000000000000e+02];
in_data.mater(9).Note='mass spring model';




end