function [S,maxF,F] = snst_maxF_dynamic_model(randInp,inp,mat_trk,mat_ws,geo,N)
parfor nsim = 1:N
    inpOne = inp;
    for npara = 1:8
        nmat = setnumC1(npara);
        ndata = setnumC2(npara);
        inpOne.mater(nmat).Data(ndata)= randInp(npara,nsim);
    end
    inpOne.mater(1).Data(5)=inpOne.mater(1).Data(1)/(2*(1+0.3));
    inpOne.mater(2).Data(5)=inpOne.mater(2).Data(1)/(2*(1+0.17));
    
    mat_trk = form_mat_trk_2(inpOne,geo);
    
    F = dynamic_main_snst_mcs(inpOne,mat_trk,geo,mat_ws);
    
    S(nsim).F=F;
    %S(npara).dis=[dis_out,dis.w];
    %S(npara).vel=[vel_out,vel.w];
    %S(npara).acc=[acc_out,acc.w];
    
    disp (['Simulation No. ',num2str(nsim),' finished. Time: ' datestr(now)]);
end