% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Sensitivity analysis using hammer test
% % Author: Chen Shen
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[inp,NNslpf] = get_input_vtrack();
% para = linspace(0.1,1.1,10);%logspace(-1,1,10);%linspace(0.1,1.1,10);%para = [0.9,1,1.1];linspace(0.1,2,20);%[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8];
% para = round(para,1);
para = [0.5:0.5:4]*1e8;
% para = [5e7:5e7:30e7]/NNslpf/5;
% para = [1:1:7]*1e4/25;
% para = [0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3];
nmat = 3;
ndata = 1; 
refval=inp.mater(nmat).Data(ndata);
btypr=inp.mesh.btypr;
btyps=inp.mesh.btyps;
[nodeCoord] = node_coor(inp,1);
S=struct('H1',zeros(length(para),1));
% S1=struct('frq',zeros(10,1));
% SfrqEigb = zeros(length(mat_trk.K_reduced),N);
% Swavenum = zeros(length(mat_trk.K_reduced),N); 
% w=genexpwin(1000);

for i = 1:length(para)
   inpOne = inp;
   inpOne.mater(nmat).Data(ndata)= para(i);
[geo] = mesh_trk_full(btypr,btyps,nodeCoord);
mat_trk = form_mat_trk_2(inpOne,geo);
% !! added eigen-analysis 
Kr=round(full(mat_trk.K_reduced),-1);
Mr=round(full(mat_trk.M_reduced),12);
if issymmetric(Kr) && issymmetric(Mr)
    [V,D]=eig(Kr,Mr);
frqEigb=abs(sqrt(diag(D)));
frqEigb=frqEigb./(2*pi());
indNodeRail = geo.ND(:,5) == 1;
x=geo.ND(indNodeRail,2);
[wavenum,fv,Vq,xq]=shape2dispersion(V(1:2:sum(indNodeRail == 1)*2,:),x',2048);

SfrqEigb(:,i) = frqEigb;
Swavenum(:,i) = wavenum(:,2);
else
    V = zeros(size(full(mat_trk.K_reduced)));
    D = eye(size(full(mat_trk.K_reduced)));
    disp(['Kr or Mr is asymmetric for simulation No. ', num2str(i)]);
end
% !!
H1 = hammer_main_snst_mcs(inpOne,mat_trk,geo);

% [H1,H2,~,pxx,pxy,fxx]=tran_fun([force(1:1000)',dis_x_load(1:1000)],w,0,25600,25000);
S(i).H1=H1;

end

% for i = 1:length(para)
% frqout=clean_up(S(i).frq(2:15),30);
% frqout=[S(i).frq(1);frqout];
% S1(i).frq=frqout;
% end
% 
% function seq_out=clean_up(seq,tol)
% diff=seq(2:end)-seq(1:end-1);
% ind=diff<tol;
% ind=[0; ind];
% seq_out=seq(~ind);
% end