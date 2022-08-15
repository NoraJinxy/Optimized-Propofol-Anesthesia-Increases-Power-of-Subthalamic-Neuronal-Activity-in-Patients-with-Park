clear all
close all
clc

load Fig_SpectrumData
%%
figure
ampli_LA_sdi=interp1(freq,ampli_LA-sd_LA,freq,'spline');
ampli_LA_sdj=interp1(freq,ampli_LA+sd_LA,freq,'spline');
ampli_GA_sdi=interp1(freq,ampli_GA-sd_GA,freq,'spline');
ampli_GA_sdj=interp1(freq,ampli_GA+sd_GA,freq,'spline');
hold on
h8=fill([freq,fliplr(freq)],[ampli_LA_sdi,fliplr(ampli_LA_sdj)],'r');
set(h8,'edgealpha',0,'facealpha',0.3) 
h9=fill([freq,fliplr(freq)],[ampli_GA_sdi,fliplr(ampli_GA_sdj)],'b');
set(h9,'edgealpha',0,'facealpha',0.3) 
plot(freq,ampli_LA,'r','linewidth',2)
plot(freq,ampli_GA,'b','linewidth',2)
y=[-10,6];
plot([2.9999,3],y,'k--');text(0.5,max(y)-1,'delta')
plot([7.9999,8],y,'k--');text(4.5,max(y)-1,'theta')
plot([13.9999,14],y,'k--');text(10.5,max(y)-1,'alpha')
text(18,max(y)-1,'beta')
g5=legend('LA','GA','mean LA','mean GA');
set(g5,'fontsize',11);
xlabel('Frequency (Hz)','fontsize',12);ylabel('Power (dB)','fontsize',12);
% title('局麻组、全麻组下丘脑核')
hold off
ylim([-10 6])
xlim([0 25])
%%

set(gcf,'PaperPositionMode','auto')
%print('-depsc2',figure(1),'Fig_Spectrum2.eps')
print('-dtiff','-r300','Fig_Spectrum')