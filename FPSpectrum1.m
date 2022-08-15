clear all
close all
clc
%% Load Data
dirname = 'FPResults\';
w = what(dirname);%列出文件下内容
fpresultsfilenames = w.mat;

%% extract index for matrix concatenation
Num_Awake = 0; Num_TCI = 0;
str_subject_Awake = {}; str_subject_TCI = {};

for n_file=1:numel(fpresultsfilenames) % START
    n_file
    filename = [dirname,fpresultsfilenames{n_file}];
    load (filename)
    Fq_All(n_file,:) = 10*log10((Fq.^2)/length(Fq)); %plot power specturm
%     Fq_All(n_file,:) = Fq; %plot frequency specturm
%     Fq_All(n_file,:)=10*log(Fq); %plot violin
    Pw_All(n_file,:)=Pw;
    f_All(n_file,:)=f;

    IDX_Group_All{n_file,1} = fpresultsfilenames{n_file}(1:2);
    IDX_Subject_All{n_file,1} = fpresultsfilenames{n_file}(4:5);
    IDX_Lateral_All{n_file,1} = fpresultsfilenames{n_file}(7);
    IDX_Depth_All{n_file,1} = fpresultsfilenames{n_file}(9:12);
end

idx_group_Awake = find(strcmp(IDX_Group_All,'LA'));
idx_group_TCI = find(strcmp(IDX_Group_All,'GA'));
idx_lateral_L = find(strcmp(IDX_Lateral_All,'L'));
idx_lateral_R = find(strcmp(IDX_Lateral_All,'R'));

str_subject_Awake = unique(IDX_Subject_All(idx_group_Awake));
Num_Awake = numel(str_subject_Awake);
str_subject_TCI = unique(IDX_Subject_All(idx_group_TCI));
Num_TCI = numel(str_subject_TCI);

% index of subject
for n_subj_Awake = 1:Num_Awake
    idx_subj_Awake{n_subj_Awake,1}=find(strcmp(IDX_Group_All,'LA')&strcmp(IDX_Subject_All,str_subject_Awake(n_subj_Awake)));
end
for n_subj_TCI = 1:Num_TCI
    idx_subj_TCI{n_subj_TCI,1}=find(strcmp(IDX_Group_All,'GA')&strcmp(IDX_Subject_All,str_subject_TCI(n_subj_TCI)));
end

%% extract data
% GA
for n_subj_TCI = 1:Num_TCI
    tmp_idx = intersect(idx_subj_TCI{n_subj_TCI,1},idx_lateral_L);
    Fq_TCI_L(n_subj_TCI,:) = mean(Fq_All(tmp_idx,:));
    Pw_TCI_L(n_subj_TCI,:) = mean(Pw_All(tmp_idx,:));
    f_TCI_L(n_subj_TCI,:) = mean(f_All(tmp_idx,:));
    tmp_idx = intersect(idx_subj_TCI{n_subj_TCI,1},idx_lateral_R);
    Fq_TCI_R(n_subj_TCI,:) = mean(Fq_All(tmp_idx,:));
    Pw_TCI_R(n_subj_TCI,:) = mean(Pw_All(tmp_idx,:));
    f_TCI_R(n_subj_TCI,:) = mean(f_All(tmp_idx,:));
end
for n_subj_Awake = 1:Num_Awake
    tmp_idx = intersect(idx_subj_Awake{n_subj_Awake,1},idx_lateral_L);
    Fq_Awake_L(n_subj_Awake,:) = mean(Fq_All(tmp_idx,:));
    Pw_Awake_L(n_subj_Awake,:) = mean(Pw_All(tmp_idx,:));
    f_Awake_L(n_subj_Awake,:) = mean(f_All(tmp_idx,:));
    tmp_idx = intersect(idx_subj_Awake{n_subj_Awake,1},idx_lateral_R);
    Fq_Awake_R(n_subj_Awake,:) = mean(Fq_All(tmp_idx,:));
    Pw_Awake_R(n_subj_Awake,:) = mean(Pw_All(tmp_idx,:));
    f_Awake_R(n_subj_Awake,:) = mean(f_All(tmp_idx,:));
end
Fq_Awake = (Fq_Awake_L+Fq_Awake_R)/2;
Fq_TCI = (Fq_TCI_L+Fq_TCI_R)/2;
Pw_Awake = (Pw_Awake_L+Pw_Awake_R)/2;
Pw_TCI = (Pw_TCI_L+Pw_TCI_R)/2;
f_Awake = (f_Awake_L+f_Awake_R)/2;
f_TCI = (f_TCI_L+f_TCI_R)/2;

%% violin plot
order1={'delta L','delta R','theta L','theta R','alpha L','alpha R','beta L','beta R','gamma L','gamma R'};
order2={'delta Awake','delta TCI','theta Awake','theta TCI','alpha Awake','alpha TCI','beta Awake','beta TCI','gamma Awake','gamma TCI'};
ColorDesign=[0,0.5,0.7;1,0,0];
[TCI_Average,TCI_Name]=Violin_Matrix(1,Num_TCI,Num_TCI,Fq_TCI_L,Fq_TCI_R,f);

[LA_Average,LA_Name]=Violin_Matrix(1,Num_Awake,Num_Awake,Fq_Awake_L,Fq_Awake_R,f);
figure
hold on
LA_vs = violinplot(LA_Average, LA_Name,ColorDesign,'GroupOrder',order1);
rectangle('Position',[6.4,7,0.5,1],'facecolor',[1,0.7,0.8]);
rectangle('Position',[6.4,5,0.5,1],'facecolor',[0.4,0.7,0.9]);
rectangle('Position',[6.3,2.5,4,6]);
plot(6.65,3.5,'r*')
text(7,7.5,'Left');
text(7,5.5,'Right');
text(7,3.5,'significant different');
ylabel('Average Frequency/dB');
title('Awake_L vs Awake_R')
hold off
xlim([0 11])
%{
figure
hold on
TCI_vs = violinplot(TCI_Average, TCI_Name,ColorDesign,'GroupOrder',order1);
rectangle('Position',[6.4,7,0.5,1],'facecolor',[1,0.7,0.8]);
rectangle('Position',[6.4,5,0.5,1],'facecolor',[0.4,0.7,0.9]);
rectangle('Position',[6.3,2.5,4,6]);
plot(6.65,3.5,'r*')
text(7,7.5,'Left');
text(7,5.5,'Right');
text(7,3.5,'significant different');
ylabel('Average Frequency/dB');
title('TCI_L vs TCI_R')
hold off
xlim([0 11])

[L_Average,L_Name]=Violin_Matrix(2,Num_TCI,Num_Awake,Fq_TCI_L,Fq_Awake_L,f);
figure
hold on
L_vs = violinplot(L_Average, L_Name,ColorDesign,'GroupOrder',order2);
rectangle('Position',[6.4,7,0.5,1],'facecolor',[1,0.7,0.8]);
rectangle('Position',[6.4,5,0.5,1],'facecolor',[0.4,0.7,0.9]);
rectangle('Position',[6.3,2.5,4,6]);
plot(6.65,3.5,'r*')
rectangle('Position',[3,5,1,0.1],'facecolor','k');
plot(3.5,6,'r*');
text(7,7.5,'Awake');
text(7,5.5,'TCI');
text(7,3.5,'significant different');
ylabel('Average Frequency/dB');
title('Awake-L vs TCI-L')
hold off
xlim([0 11])

[R_Average,R_Name]=Violin_Matrix(2,Num_TCI,Num_Awake,Fq_TCI_R,Fq_Awake_R,f);
figure
hold on
R_vs = violinplot(R_Average, R_Name,ColorDesign,'GroupOrder',order2);
rectangle('Position',[6.4,7,0.5,1],'facecolor',[1,0.7,0.8]);
rectangle('Position',[6.4,5,0.5,1],'facecolor',[0.4,0.7,0.9]);
rectangle('Position',[6.3,2.5,4,6]);
plot(6.65,3.5,'r*')
rectangle('Position',[1,7,1,0.1],'facecolor','k');
plot(1.5,8,'r*');
text(7,7.5,'Awake');
text(7,5.5,'TCI');
text(7,3.5,'significant different');
ylabel('Average Frequency/dB');
title('Awake-R vs TCI-R')
hold off
xlim([0 11])

[TCIAwake_Average,LR_Name]=Violin_Matrix(2,Num_TCI,Num_Awake,Fq_TCI,Fq_Awake,f);
figure
hold on
LR_vs = violinplot(TCIAwake_Average, LR_Name,ColorDesign,'GroupOrder',order2);
rectangle('Position',[6.4,7,0.5,1],'facecolor',[1,0.7,0.8]);
rectangle('Position',[6.4,5,0.5,1],'facecolor',[0.4,0.7,0.9]);
rectangle('Position',[6.3,2.5,4,6]);
plot(6.65,3.5,'r*')
rectangle('Position',[1,4.5,1,0.1],'facecolor','k');
rectangle('Position',[3,3.5,1,0.1],'facecolor','k');
plot(1.5,5.5,'r*');
plot(3.5,4.5,'r*');
text(7,7.5,'Awake');
text(7,5.5,'TCI');
text(7,3.5,'significant different');
ylabel('Average Frequency/dB');
title('TCI vs Awake');
hold off
xlim([0 11])
%}
%% frequency specturm
db_Fq_Awake_L=Fq_Awake_L(:,6:end);db_Fq_Awake_R=Fq_Awake_R(:,6:end);
db_Fq_TCI_L=Fq_TCI_L(:,6:end);db_Fq_TCI_R=Fq_TCI_R(:,6:end);
db_Fq_Awake=Fq_Awake(:,6:end);db_Fq_TCI=Fq_TCI(:,6:end);

% db_Fq_LA_L=Fq_LA_L;db_Fq_LA_R=Fq_LA_R;
% db_Fq_GA_L=Fq_GA_L;db_Fq_GA_R=Fq_GA_R;
% db_Fq_LA=Fq_LA;db_Fq_GA=Fq_GA;

m_Fq_Awake_L = mean(db_Fq_Awake_L);
sd_Fq_Awake_L = std(db_Fq_Awake_L);
m_Fq_Awake_R = mean(db_Fq_Awake_R);
sd_Fq_Awake_R = std(db_Fq_Awake_R);
m_Fq_TCI_L = mean(db_Fq_TCI_L);
sd_Fq_TCI_L = std(db_Fq_TCI_L);
m_Fq_TCI_R = mean(db_Fq_TCI_R);
sd_Fq_TCI_R = std(db_Fq_TCI_R);
freq=f(6:end);
% freq=f;
%{
figure
ampli_Awake_L_sdi=interp1(freq,m_Fq_Awake_L-sd_Fq_Awake_L,freq,'spline');
ampli_Awake_L_sdj=interp1(freq,m_Fq_Awake_L+sd_Fq_Awake_L,freq,'spline');
ampli_Awake_R_sdi=interp1(freq,m_Fq_Awake_R-sd_Fq_Awake_R,freq,'spline');
ampli_Awake_R_sdj=interp1(freq,m_Fq_Awake_R+sd_Fq_Awake_R,freq,'spline');
hold on
h=fill([freq,fliplr(freq)],[ampli_Awake_L_sdi,fliplr(ampli_Awake_L_sdj)],'r')
set(h,'edgealpha',0,'facealpha',0.3) 
h1=fill([freq,fliplr(freq)],[ampli_Awake_R_sdi,fliplr(ampli_Awake_R_sdj)],'b')
set(h1,'edgealpha',0,'facealpha',0.3) 
plot(freq,m_Fq_Awake_L,'r','linewidth',2)
plot(freq,m_Fq_Awake_R,'b','linewidth',2)
y=[-12,6];
plot([2.9999,3],y,'k');text(0.5,min(y)+1,'delta')
plot([7.9999,8],y,'k');text(4.5,min(y)+1,'theta')
plot([13.9999,14],y,'k');text(10.5,min(y)+1,'alpha')
text(18,min(y)+1,'beta')
g1=legend('Awake_L SD','Awake_R SD','mean Awake_L','mean Awake_R');
set(g1,'fontsize',11);
xlabel('Frequency(Hz)','fontsize',18);ylabel('Power','fontsize',18);
% title('局麻组左右下丘脑核')
hold off
% ylim([10,70])
xlim([0 25])

figure
ampli_TCI_L_sdi=interp1(freq,m_Fq_TCI_L-sd_Fq_TCI_L,freq,'spline');
ampli_TCI_L_sdj=interp1(freq,m_Fq_TCI_L+sd_Fq_TCI_L,freq,'spline');
ampli_TCI_R_sdi=interp1(freq,m_Fq_TCI_R-sd_Fq_TCI_R,freq,'spline');
ampli_TCI_R_sdj=interp1(freq,m_Fq_TCI_R+sd_Fq_TCI_R,freq,'spline');
hold on
h2=fill([freq,fliplr(freq)],[ampli_TCI_L_sdi,fliplr(ampli_TCI_L_sdj)],'r')
set(h2,'edgealpha',0,'facealpha',0.3) 
h3=fill([freq,fliplr(freq)],[ampli_TCI_R_sdi,fliplr(ampli_TCI_R_sdj)],'b')
set(h3,'edgealpha',0,'facealpha',0.3) 
plot(freq,m_Fq_TCI_L,'r','linewidth',2)
plot(freq,m_Fq_TCI_R,'b','linewidth',2)
y=[-12,6];
plot([2.9999,3],y,'k');text(0.5,min(y)+1,'delta')
plot([7.9999,8],y,'k');text(4.5,min(y)+1,'theta')
plot([13.9999,14],y,'k');text(10.5,min(y)+1,'alpha')
text(18,min(y)+1,'beta')
g2=legend('TCI_L SD','TCI_R SD','mean TCI_L','mean TCI_R');
set(g2,'fontsize',11);
xlabel('Frequency(Hz)','fontsize',18);ylabel('Power','fontsize',18);
% title('全麻组左右下丘脑核')
hold off
% ylim([-5,15])
xlim([0 25])

figure
ampli_TCI_L_sdi=interp1(freq,m_Fq_TCI_L-sd_Fq_TCI_L,freq,'spline');
ampli_TCI_L_sdj=interp1(freq,m_Fq_TCI_L+sd_Fq_TCI_L,freq,'spline');
ampli_Awake_L_sdi=interp1(freq,m_Fq_Awake_L-sd_Fq_Awake_L,freq,'spline');
ampli_Awake_L_sdj=interp1(freq,m_Fq_Awake_L+sd_Fq_Awake_L,freq,'spline');
hold on
h4=fill([freq,fliplr(freq)],[ampli_Awake_L_sdi,fliplr(ampli_Awake_L_sdj)],'r')
set(h4,'edgealpha',0,'facealpha',0.3) 
h5=fill([freq,fliplr(freq)],[ampli_TCI_L_sdi,fliplr(ampli_TCI_L_sdj)],'b')
set(h5,'edgealpha',0,'facealpha',0.3) 
plot(freq,m_Fq_Awake_L,'r','linewidth',2)
plot(freq,m_Fq_TCI_L,'b','linewidth',2)
y=[-12,6];
plot([2.9999,3],y,'k');text(0.5,min(y)+1,'delta')
plot([7.9999,8],y,'k');text(4.5,min(y)+1,'theta')
plot([13.9999,14],y,'k');text(10.5,min(y)+1,'alpha')
text(18,min(y)+1,'beta')
g3=legend('Awake_L SD','TCI_L SD','mean Awake_L','mean TCI_L');
set(g3,'fontsize',11);
xlabel('Frequency(Hz)','fontsize',18);ylabel('Power','fontsize',18);
% title('局麻、全麻组左下丘脑核')
hold off
% ylim([-5,15])
xlim([0 25])

figure
ampli_TCI_R_sdi=interp1(freq,m_Fq_TCI_R-sd_Fq_TCI_R,freq,'spline');
ampli_TCI_R_sdj=interp1(freq,m_Fq_TCI_R+sd_Fq_TCI_R,freq,'spline');
ampli_Awake_R_sdi=interp1(freq,m_Fq_Awake_R-sd_Fq_Awake_R,freq,'spline');
ampli_Awake_R_sdj=interp1(freq,m_Fq_Awake_R+sd_Fq_Awake_R,freq,'spline');
hold on
h6=fill([freq,fliplr(freq)],[ampli_Awake_R_sdi,fliplr(ampli_Awake_R_sdj)],'r')
set(h6,'edgealpha',0,'facealpha',0.3) 
h7=fill([freq,fliplr(freq)],[ampli_TCI_R_sdi,fliplr(ampli_TCI_R_sdj)],'b')
set(h7,'edgealpha',0,'facealpha',0.3) 
plot(freq,m_Fq_Awake_R,'r','linewidth',2)
plot(freq,m_Fq_TCI_R,'b','linewidth',2)
y=[-12,6];
plot([2.9999,3],y,'k');text(0.5,min(y)+1,'delta')
plot([7.9999,8],y,'k');text(4.5,min(y)+1,'theta')
plot([13.9999,14],y,'k');text(10.5,min(y)+1,'alpha')
text(18,min(y)+1,'beta')
g4=legend('Awake_R SD','TCI_R SD','mean Awake_R','mean TCI_R');
set(g4,'fontsize',11);
xlabel('Frequency(Hz)','fontsize',18);ylabel('Power','fontsize',18);
% title('局麻、全麻组右下丘脑核')
hold off
% ylim([-2,10])
xlim([0 25])
%}
%% frequency spectrum TCI vs Awake
ampli_Awake=mean(db_Fq_Awake);
ampli_TCI=mean(db_Fq_TCI);
sd_Awake=std(db_Fq_Awake);
sd_TCI=std(db_Fq_TCI);

figure
ampli_Awake_sdi=interp1(freq,ampli_Awake-sd_Awake,freq,'spline');
ampli_Awake_sdj=interp1(freq,ampli_Awake+sd_Awake,freq,'spline');
ampli_TCI_sdi=interp1(freq,ampli_TCI-sd_TCI,freq,'spline');
ampli_TCI_sdj=interp1(freq,ampli_TCI+sd_TCI,freq,'spline');
hold on
h8=fill([freq,fliplr(freq)],[ampli_Awake_sdi,fliplr(ampli_Awake_sdj)],'r')
set(h8,'edgealpha',0,'facealpha',0.3) 
h9=fill([freq,fliplr(freq)],[ampli_TCI_sdi,fliplr(ampli_TCI_sdj)],'b')
set(h9,'edgealpha',0,'facealpha',0.3) 
plot(freq,ampli_Awake,'r','linewidth',2)
plot(freq,ampli_TCI,'b','linewidth',2)
y=[-12,6];
plot([2.9999,3],y,'k');text(0.5,min(y)+1,'delta')
plot([7.9999,8],y,'k');text(4.5,min(y)+1,'theta')
plot([13.9999,14],y,'k');text(10.5,min(y)+1,'alpha')
text(18,min(y)+1,'beta')
g5=legend('Awake SD','TCI SD','mean Awake','mean TCI');
set(g5);
xlabel('Frequency(Hz)');ylabel('Power');
% title('局麻组、全麻组下丘脑核')
hold off
% ylim([-2 10])
% xlim([0 25])
