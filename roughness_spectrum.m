%%%%%%%
% Generate rail roughness spectrum based on 6th grade US (??)
% Author: Pan Zhang
% 
%%%%
function [Xnz,t] = roughness_spectrum(Vmax, dt, Nr)
% clear all
% clc
% clf
%% �����ƽ˳����
% Vmax=13/3.6;                                              %�г���������ٶȣ�m/s)
BCmin=0.005;                                               %��ƽ˳��С������m��
BCmax=100;                                              %��ƽ˳��󲨳���m��
fmin=Vmax/BCmax;                                         %��СƵ�ʣ�Hz��
fmax=Vmax/BCmin;                                         %���Ƶ�ʣ�Hz��

%% ��������
% dt=1/25000;                                                %��������(s)
% Nr=4400;%2^16;                                                %ʱ���������
Ts=dt*Nr;                                                   %ģ����ʱ�䣨s��
df=1/(Nr*dt);                                            %%
Nf=round((fmax-fmin)/df);                                %��ЧƵ�ʶεĲ�������
No=round(fmin/df);                                       %%

%������������׸ߵͲ�ƽ˳����
% k=0.25;                                                  %��ȫϵ��
% Av=0.0339*Vmax/(2*pi);                                               %�ֲڶȳ�����cm2*Hz)
% Oc=0.8245*Vmax/(2*pi);                                               %�ض�Ƶ�ʣ�Hz��

%������������׸ߵͲ�ƽ˳����
k=0.25;                                                  %��ȫϵ��
Av=0.2095*Vmax/(2*pi);                                               %�ֲڶȳ�����cm2*Hz)
Oc=0.8245*Vmax/(2*pi);                                               %�ض�Ƶ�ʣ�Hz��
%Av=0.0339*2*pi/Vmax;                         %�ֲڶȳ�����cm2*Hz)
%�������ƽ˳ʱ��ģ������
for m=No:1:Nf
f=fmin+(m-No)*df;  
% f=m*df;
fb(m-No+1)=f;
Sv(m-No+1)=(k.*Av.*Oc.*Oc)./(f.*f.*f.*f+f.*f.*Oc.*Oc)/1e4;           %���������׹���ߵͲ�ƽ˳
% if f<Vmax
% Sv(m-No+1)=(k.*Av.*Oc.*Oc)./(f.^4+f.^2*Oc.*Oc);            %���������׹���ߵͲ�ƽ˳
% else
% F=f/Vmax;   
% Sv(m-No+1)=0.036*F.^(-3.15)/(100*Vmax);
% end
end
% m=No:1:Nf;
% f=m*df;
% Sv=(k.*Av.*Oc.*Oc)./(f.*f.*f.*f+f.*f.*Oc.*Oc);           %���������׹���ߵͲ�ƽ˳
absXk=Nr*sqrt(Sv*df/2);                                    %ʱ�����е�Ƶ��ģֵ
J=2*pi*rand(size(Sv));                                   %����0-2PI�ľ��ȷֲ�
E=exp(1i*J);                                              %������λ����
Xk1=zeros(1,No);                                         %��������������ĿΪNo
Xk2=E.*absXk;                                            %��ЧƵ�ʶε�Ƶ��
Xk3=zeros(1,(Nr/2-Nf));                                  %��������������ĿΪ(Nr/2-Nf+1)
Xk4=[Xk1,Xk2,Xk3];                                       %(0-Nr/2)���е�Ƶ��
Xk41=conj(Xk4);                                          %ȡ(0-Nr/2)���е�Ƶ�׵Ĺ���
Xk5=fliplr(Xk41);                                         %��Ƶ��Xk4���ڵ�Nr/2���ԳƱ任
Xk5(1)=[];
Xk5(Nr/2)=[];
Xk=[Xk4,Xk5];                                            %(0,Nr)���е�Ƶ��
xn=ifft(Xk);                                             %����Ҷ��任
Xn=real(xn);                                             %ȡʵ������任���鲿��Ϊ0
Xnz=Xn.';                                                %��ʱ������ת��Ϊ������
n=0:1:(Nr-1);                                            %%
t=n*dt;                                                  %����ʱ������
tz=t.'*Vmax;                                                  %��ʱ������ת��Ϊ������
Xnt=[tz,Xnz];                                            %���ʱ�����к�ʱ������ 
%f=Vmax./f;
ft=f.';                                                  %��ʱ������ת��Ϊ������
%Sv=2*pi*Sv/(Vmax);
Svt=Sv.';
% Svf=[ft,Svt];                                            %���ʱ�����к�ʱ������
% save('Svf.txt','Svf','-ascii');                          %�����ݱ������ı�Xnt��
%save('guidao2.txt','Svt','-ascii');                          %�����ݱ������ı�Xnt��

%��������
% save('����5����160.txt','Xnt','-ascii');                          %�����ݱ������ı�Xnt��
% figure;plot(t,Xnz)
