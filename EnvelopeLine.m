clear all
% close all
clc
tic
%% LA
% load LA_14_L_-02500_0017
% fs = 24000; % sampling rate
% s= x.signal; % signal
% t = x.time/fs; % time
% T = max(t); % data length
% h=hilbert(s);

% time=t(1:fs);
% y=s(9*fs+1:10*fs);
% h1=h(9*fs+1:10*fs);

%% GA
load GA_02_L_-03000_0018
% load GA_15_R_+00000_0043
s= x.signal;
fs = 24000; % sampling rate
t = x.time/fs; % time
T = max(t); % data length
h=hilbert(s);
y=s(1:fs);
h1=h(1:fs);
time=t(1:fs);


%% spike selected

noise_std = median(abs(s)/0.6745);
TH = 4*noise_std;
spk_idx = find((y)>TH);
if isempty(spk_idx)
spike_train = zeros(size(time));
MBI = NaN;
FR = NaN;
else
d_spk_idx = diff(spk_idx);
idx_seg = find(d_spk_idx>1);%find---返回位置

N_Spk = numel(idx_seg)+1;
peak_idx = [];
for n=1:N_Spk
    if n==1
        spk_idx_c{n,1}=spk_idx([1:idx_seg(n)]);
    elseif n<N_Spk
        spk_idx_c{n,1}=spk_idx([(idx_seg(n-1)+1):idx_seg(n)]);
    else
        spk_idx_c{n,1}=spk_idx([(idx_seg(n-1)+1):numel(spk_idx)]);
    end
    [tmp,peak_idx_tmp] = max(y(spk_idx_c{n,1}));
    peak_idx(n,1)=spk_idx_c{n,1}(peak_idx_tmp);
end
spike_train = zeros(size(t));
spike_train(peak_idx) = 1;
ISI = diff(peak_idx)/fs; % ISI
MBI = numel(find(ISI<=0.01))/numel(find(ISI)>0.01); % 10ms = 0.01s
FR = numel(ISI)/T;
end % isempty


% figure
% % subplot(231)
% hold on 
% plot(time,y,'k')
% plot(time(peak_idx),y(peak_idx),'r.')
% plot(0.35,200,'r*','linewidth',5)
% set(gca,'ytick',-200:100:300)
% set(gca,'xtick',0:0.2:1)
% xlabel('Time(s)');
% ylabel('Amplitude(μV)')
% ylim([-200,300])

figure
% subplot(232)
hold on 
plot(time,y,'k')
plot(time(peak_idx),y(peak_idx),'r.')
plot(time,abs(h1(1:24000)),'r')
plot(0.35,200,'r*','linewidth',5)
xlim([0.34,0.36])
xlabel('Time(s)');
ylabel('Amplitude(μV)')
ylim([-200,300])
set(gca,'ytick',-200:100:300)


%% Autocorrelogram
% window_width=0.25*fs; %150ms
% slide=window_width/2;
% Begin=1;
% End=window_width;
% slide_num=(T*fs-window_width)/slide;
% signal=s(Begin:End); %取长度为150ms的信号
% N=numel(signal);
% r=numel(peak_idx)/(T*fs); 
% for i=1:slide_num+1
%     signal=s(Begin+slide*(i-1):End+slide*(i-1));
%     for m=1:N-1
%         c_j=[];
%         for n=1:N-m
%             c_j(n)=(signal(n)-r)*(signal(n+m)-r);
%         end
%         n_c=numel(c_j);
%         y_r(m)=sum(c_j)/n_c;
%     end
%     correlation(i,:)=y_r;
% end
% % center=sum(signal.*signal)/24000;
% Autoccorelation=mean(correlation);
% center=0;
% x_r=1-N:N-1;
% x_r=x_r/24;
% Autoccorelation=[fliplr(Autoccorelation),center,Autoccorelation];
% % y_r=y_r/center;
% figure
% plot(x_r,Autoccorelation,'k');%xlim([-150 150])
% xlabel('Time/ms');
% ylim([-200 1000])
% set(gca,'ytick',-200:200:1000)
% box off

signal=s;
N=numel(signal);
r=numel(peak_idx)/(T*fs); 
for m=1:N-1
    c_j=[];
    for n=1:N-m
        c_j(n)=(signal(n)-r)*(signal(n+m)-r);
    end
    n_c=numel(c_j);
    y_r(m)=sum(c_j)/n_c;
end
Autoccorelation=y_r;
center=0;
x_r=1-N:N-1;
x_r=x_r/24;
Autoccorelation=[fliplr(Autoccorelation),center,Autoccorelation];
% y_r=y_r/center;
figure
plot(x_r,Autoccorelation,'k');%xlim([-150 150])
xlabel('Time/ms');
ylim([-200 1000])
set(gca,'ytick',-200:200:1000)
box off
toc