clear all

%% Load Data
dirname = 'MATFilesSelected\';
w = what(dirname);
matfilenames = w.mat;

for n_file=1:numel(matfilenames) % START
    
filename = [dirname,matfilenames{n_file}];
load (filename)

disp(matfilenames{n_file})

str_group = matfilenames{n_file}(1:2);
str_subject = matfilenames{n_file}(4:5);
str_lateral = matfilenames{n_file}(7);
str_depth = matfilenames{n_file}(9:12);

s = x.signal; % signal
fs = 24000; % sampling rate
t = x.time/fs; % time  
T = max(t); % data length

%% Filtering (NOT NECESSARY)
% f_bp = [500,1000];
% n_order = 2; Wn = f_bp/(fs/2);
% [b,a] = butter(n_order,Wn);
% y = filtfilt(b,a,s);
y = s; % no filter

%% Spike Detection and Parameter Extraction
noise_std = median(abs(y)/0.6745);
TH = 4*noise_std;
spk_idx = find((y)>TH);
if isempty(spk_idx)
spike_train = zeros(size(t));
MBI = NaN;
FR = NaN;
else
d_spk_idx = diff(spk_idx);
idx_seg = find(d_spk_idx>1);
ISI = d_spk_idx(idx_seg)/fs; % ISI
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
MBI = numel(find(ISI<=0.01))/numel(find(ISI)>0.01); % 10ms = 0.01s
FR = numel(ISI)/T;
end % isempty
%{
%% Autocorrelogram
window_width=12000; %30s
slide=window_width/2;
Begin=1;
End=12000;
slide_num=(T*fs-window_width)/slide;
signal=s(Begin:End); %取长度为1s的信号
N=numel(signal);
r=numel(peak_idx)/(T*fs); 
for i=1:slide_num+1
    signal=s(Begin+slide*(i-1):End+slide*(i-1));
    for m=1:N-1
        c_j=[];
        for n=1:N-m
            c_j(n)=(signal(n)-r)*(signal(n+m)-r);
        end
        n_c=numel(c_j);
        y_r(m)=sum(c_j)/n_c;
    end
    correlation(i,:)=y_r;
end
% center=sum(signal.*signal)/24000;
Autoccorelation=mean(correlation);
center=0;
x_r=1-N:N-1;
x_r=x_r/24;
Autoccorelation=[fliplr(Autoccorelation),center,Autoccorelation];
% y_r=y_r/center;
figure
set(gcf,'Visible','off');
plot(x_r,Autoccorelation);xlim([-150 150]);
title('autocorrelogram','fontsize',14);
xlabel('Time/ms','fontsize',14);

figurename1 = ['Autocorrelograms\',matfilenames{n_file}(1:end-4)];
print('-dtiff','-r200',figurename1)
close all
%}
%% Show Data
figure
% set(gcf,'Visible','off');
subplot(311)
t=1:24000;
s=s(1:24000);
plot(t,s,'k')
xlabel('Time (s)'); ylabel('Amplitude (\it\muV)');
title([str_group,' | ',str_subject,' | ',str_lateral,' | ',str_depth],'fontsize',12);
axis tight
subplot(312)
hold on; box on;
plot(t,y,'k')
plot(t(peak_idx),y(peak_idx),'r.')
xlabel('Time (s)'); ylabel('Amplitude (\it\muV)');
axis tight
subplot(313)
stem(t,spike_train,'r.')
xlabel('Time (s)'); ylabel('Amplitude (\it\muV)');
axis tight

% figurename = ['Figures\',matfilenames{n_file}(1:end-4)];
% print('-dtiff','-r200',figurename)
% close all


%% save file
resultname = ['Results\',matfilenames{n_file}(1:end-4)];
save (resultname,'spike_train','MBI','FR')

%% spike_train frequency and power spectrum
%Frequency spectrum

signal=spike_train;
L=length(signal);
f=fs*(0:(L/2))/L;
f=f(1:1001);
ft=fft(signal,L);
fq=abs(ft);
Fq=fq(1:L/2+1);
Fq=Fq(1:1001);
figure
subplot(311)
plot(f,Fq);
% xlim([0 1000]);
title('Frequency spectrum');
xlabel('f(Hz)');ylabel('|Amplitude|');

%Power spectrum
Pw=10*log10((Fq.^2)/length(Fq));
%Fq.^2/length(Fq);

%% save file
resultname = ['FPResults\',matfilenames{n_file}(1:end-4)];
save (resultname,'f','Fq','Pw')

end % END OF FILES