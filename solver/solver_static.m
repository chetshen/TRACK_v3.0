%static solver, nonlinear with contact iteration
function [dis,Z,F,I_ss]=solver_static(mat_trk,inp,shape,contactID, mat_vhcl)
F0=(inp.mater(6).Data(2)+inp.mater(6).Data(1))*9.8+inp.ext_force.wh_ld;
F=[F0,F0];
switch contactID
    case 5
%-----non-linear---
penetration=(F0/inp.mater(contactID).Data)^(2/3);
    case 10
%-----linear-----
penetration=F0/inp.mater(contactID).Data;
% penetration=(F/inp.mater(contactID).Data)^(2/3);
    case 1000
        delta_ss = 0.015;
        C = 0.05;
         I_ss = sqrt(F0*delta_ss.^2/C);
penetration = delta_ss;
end


R_trk=-F*shape;
dis.r=(mat_trk.K_reduced\R_trk')';
Z.r=shape*dis.r';
Z.w=Z.r-penetration;
dis.w=Z.w;


end
