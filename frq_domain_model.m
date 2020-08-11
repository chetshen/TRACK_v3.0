function [H_sys, recep_con, recep_w] = frq_domain_model(fxx,recep_trk,k_con,m_ws)
omega = 2*pi()*fxx;
recep_con = 1./k_con.*ones(length(fxx),1);
recep_w = 1./(omega.^2.*m_ws);
H_sys = 1./(recep_trk+recep_con+recep_w);

