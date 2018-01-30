function [H,p,res,meq,stfeq]=eig2frf(phi,wn,dr,wi)
H=zeros(length(phi),length(phi),length(wi));
p=zeros(length(wn),1);
res=zeros(length(phi),length(phi),length(wn));
paug=zeros(length(phi),length(phi),length(wn));
meq=zeros(length(phi),length(phi),length(wn));

for r=1:length(wn)
    for k=1:length(phi)
        for j=1:length(phi)
            res(j,k,r)=phi(j,r)*phi(k,r)/(2i*wn(r));
            meq(j,k,r)=1./(phi(j,r)*phi(k,r));
        end
    end
    p(r,1)=-1*wn(r)*dr(r)+wn(r)*sqrt(1-dr(r).^2)*1i;
    paug(:,:,r)=p(r,1)*ones(length(phi),length(phi));
end

stfeq=meq.*reshape(wn,1,1,[]).^2;

for w=1:length(wi)
    H(:,:,w)=sum(res./(wi(w)*1i-paug)+conj(res)./(wi(w)*1i-conj(paug)),3);  
end

%figure;
%frq=wi./(2*pi);
%loglog(frq,abs(squeeze(H(1,1,:))));
%
end