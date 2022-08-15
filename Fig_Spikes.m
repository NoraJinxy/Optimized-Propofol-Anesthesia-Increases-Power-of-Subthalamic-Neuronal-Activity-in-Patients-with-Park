clear all
close all
clc

load GA_02_L_-02000_0020

s = x.signal; % signal
fs = 24000; % sampling rate
t = x.time/fs; % time
T = max(t); % data length

%% spike selected
noise_std = median(abs(s)/0.6745);
TH = 4*noise_std;
spk_idx = find((s)>TH);
if isempty(spk_idx)
spike_train = zeros(size(t));
MBI = NaN;
FR = NaN;
else
d_spk_idx = diff(spk_idx);
idx_seg = find(d_spk_idx>1);%find---∑µªÿŒª÷√

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
    [tmp,peak_idx_tmp] = max(s(spk_idx_c{n,1}));
    peak_idx(n,1)=spk_idx_c{n,1}(peak_idx_tmp);
end
spike_train = zeros(size(t));
spike_train(peak_idx) = 1;
ISI = diff(peak_idx)/fs; % ISI
MBI = numel(find(ISI<=0.01))/numel(find(ISI)>0.01); % 10ms = 0.01s
FR = numel(ISI)/T;
end % isempty

%% Show Data
figure
subplot(311)
plot(t,s,'k')
xlabel('Time (s)'); 
ylabel('\it\muV','fontsize',12);
title('Raw Data','fontsize',12);
xlim([0 1])
% axis tight
subplot(312)
hold on; box on;
plot(t,s,'k')
plot(t(peak_idx),s(peak_idx),'r.')
xlabel('Time (s)'); 
ylabel('\it\muV','fontsize',12);
title('Detected Spikes','fontsize',12);
xlim([0 1])
% axis tight

subplot(313)
%stem(t,spike_train,'k--')
plot(t,spike_train,'k')
xlabel('Time (s)'); 
title('Spike Train','fontsize',12);
xlim([0 1]); ylim([0 1.5])
set(gca,'ytick',[])
% axis tight

set(gcf,'PaperPositionMode','auto')
print('-depsc2',figure(1),'Fig_spikes.eps')