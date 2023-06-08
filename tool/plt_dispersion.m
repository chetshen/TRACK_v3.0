function plt_dispersion(frqEigb_bound,frqEigb,Vfft,ind_fv,wavenum_plt)
for ii = 1:size(frqEigb_bound,1)
ind_frq = frqEigb < frqEigb_bound(ii,2) & frqEigb > frqEigb_bound(ii,1);
frq_plt = frqEigb(ind_frq);
Vfft_plt = Vfft(ind_fv,ind_frq);

wavenum_plt_twoside = [flipud(-1*wavenum_plt);wavenum_plt];
Vfft_plt_twoside = [flipud(Vfft_plt);Vfft_plt];

[wavenum_msh,frq_msh] = meshgrid(wavenum_plt_twoside,frq_plt);

mesh(wavenum_msh,frq_msh,log(Vfft_plt_twoside'))
view([0 90])
hold on
end

end