%static solver, nonlinear with contact iteration
function [dis,Z,F]=solver_static(mat_trk,inp,shape,contactID, mat_vhcl)
F=(inp.mater(6).Data(2)+inp.mater(6).Data(1))*9.8+inp.ext_force.wh_ld;
switch contactID
    case 5
%-----non-linear---
penetration=(F/inp.mater(contactID).Data)^(2/3);
    case 10
%-----linear-----
penetration=F/inp.mater(contactID).Data;
% penetration=(F/inp.mater(contactID).Data)^(2/3);
end


R_trk=-F*shape';
dis.r=(mat_trk.K_reduced\R_trk)';
Z.r=shape*dis.r';
Z.w=Z.r-penetration;
dis.w=Z.w;


end
