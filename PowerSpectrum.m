clear all
close all
clc

load PowerSpectrumData

figure
ampli_LA_sdi=interp1(freq,ampli_LA-sd_LA,freq,'spline');
ampli_LA_sdj=interp1(freq,ampli_LA+sd_LA,freq,'spline');
ampli_GA_sdi=interp1(freq,ampli_GA-sd_GA,freq,'spline');
ampli_GA_sdj=interp1(freq,ampli_GA+sd_GA,freq,'spline');
hold on
h8=fill([freq,fliplr(freq)],[ampli_LA_sdi,fliplr(ampli_LA_sdj)],'r')
set(h8,'edgealpha',0,'facealpha',0.3) 
h9=fill([freq,fliplr(freq)],[ampli_GA_sdi,fliplr(ampli_GA_sdj)],'b')
set(h9,'edgealpha',0,'facealpha',0.3) 
plot(freq,ampli_LA,'r')
plot(freq,ampli_GA,'b')
y=[-12,6];
plot([2.9999,3],y,'k');text(0.5,min(y)+1,'delta')
plot([7.9999,8],y,'k');text(4.5,min(y)+1,'theta')
plot([13.9999,14],y,'k');text(10.5,min(y)+1,'alpha')
text(18,min(y)+1,'beta')
g5=legend('LA SD','GA SD','mean LA','mean GA');
set(g5,'fontsize',11);
xlabel('Frequency(Hz)','fontsize',18);ylabel('Power','fontsize',18);
% title('局麻组、全麻组下丘脑核')
hold off
% ylim([-2 10])
xlim([0 25])