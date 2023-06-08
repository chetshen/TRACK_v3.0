function [wavenum,fv,Vq,xq,Vfft]=shape2dispersion(V,x,nq,waveL_ub,waveL_lb)
if nargin < 3
    nq=1024;
    waveL_ub = 1e9;%1.21;
    waveL_lb = 1e-9;%0.59;
elseif nargin == 3
    waveL_ub = 1e9;%1.21;
    waveL_lb = 1e-9;%0.59;
end

xq=linspace(x(1),x(end),nq);
Vq=interp1(x,V,xq);
[Vfft,fv]=pwelch(Vq,nq,0,nq,1/mean(diff(xq)));
Vfft(:,max(Vfft)<1e-20) = nan;
% [~,Ivfft]=max(Vfft);
%%
% waveL_ub = 1e9;%1.21;
% waveL_lb = 1.19;%0.59;
ind_fv = fv > 1./waveL_ub & fv < 1./waveL_lb;%155:166;
fv(~ind_fv,1) = nan;
% [~,Ivfft(index)]=max(Vfft(fv>1,index));
Vfft(~ind_fv,:) = nan;
[~,Ivfft]=max(Vfft);
%%
wavenum(:,1)=fv(Ivfft(1,:)); % wave number in m-1
%%
% wavenum(index,1)=fv(Ivfft2(1,index));
%%
wavenum(:,2)=2*pi*wavenum(:,1); % wave number in rad/m
wavenum(:,3)=1./wavenum(:,1); %wavelength in m
end

