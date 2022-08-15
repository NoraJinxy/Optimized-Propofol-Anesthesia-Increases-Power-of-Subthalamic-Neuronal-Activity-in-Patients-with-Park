clear all
close all
clc
%% Load Data
dirname = 'FPResults\';
w = what(dirname);%列出文件下内容
fpresultsfilenames = w.mat;

%% extract index for matrix concatenation
Num_LA = 0; Num_GA = 0;
str_subject_LA = {}; str_subject_GA = {};
age_LA=[64;55;61;62;44;63;64;68;55;62;42;49;66;72;57;72;61;51;55;55;64;58;...
    48;60;34;62;54;56;50;64;38;45;49];
age_GA=[53;66;52;71;36;68;52;52;65;49;50;65;43;65;69;55;43;42;58;55;47;52;...
    50;60;50;63;61;62;60;65;61;41;60;67;64];
score_pre_GA=[39;61;53;57;31;44;74;67;38;55;55;67;39;82;49;40;34;...
    46;50;54;51;46;26;71;69;50;30;48];
score_post_GA=[21;40;19;21;20;11;32;21;16;30;16;34;29;35;20;14;...
    12;14;11;18;5;14;10;26;34;16;10;20];
score_pre_LA=[41;57;54;73;61;82;68;71;97;59;93;60;58;59;...
    57;49;57;52;53;38;34;52;34;52;62;47;];
score_post_LA=[12;12;26;15;11;37;25;17;27;14;29;35;8;19;...
    22;16;15;18;22;26;22;13;20;24;16;9];
delet_LA=[13,14,15,16,17,20,31];
delet_GA=[3,4,15,24,27,29,32];

for n_file=1:numel(fpresultsfilenames) % START
    n_file
    filename = [dirname,fpresultsfilenames{n_file}];
    load (filename)
    Fq_All(n_file,:) = Fq; %frequency
%     Fq_All(n_file,:) = 10*log10((Fq.^2)/length(Fq)); %power
    Pw_All(n_file,:) = Pw;
    f_All(n_file,:)=f;

    IDX_Group_All{n_file,1} = fpresultsfilenames{n_file}(1:2);
    IDX_Subject_All{n_file,1} = fpresultsfilenames{n_file}(4:5);
    IDX_Lateral_All{n_file,1} = fpresultsfilenames{n_file}(7);
    IDX_Depth_All{n_file,1} = fpresultsfilenames{n_file}(9:12);
end

idx_group_LA = find(strcmp(IDX_Group_All,'LA'));
idx_group_GA = find(strcmp(IDX_Group_All,'GA'));
idx_lateral_L = find(strcmp(IDX_Lateral_All,'L'));
idx_lateral_R = find(strcmp(IDX_Lateral_All,'R'));

str_subject_LA = unique(IDX_Subject_All(idx_group_LA));
Num_LA = numel(str_subject_LA);
str_subject_GA = unique(IDX_Subject_All(idx_group_GA));
Num_GA = numel(str_subject_GA);

% index of subject
for n_subj_LA = 1:Num_LA
    idx_subj_LA{n_subj_LA,1}=find(strcmp(IDX_Group_All,'LA')&strcmp(IDX_Subject_All,str_subject_LA(n_subj_LA)));
end
for n_subj_GA = 1:Num_GA
    idx_subj_GA{n_subj_GA,1}=find(strcmp(IDX_Group_All,'GA')&strcmp(IDX_Subject_All,str_subject_GA(n_subj_GA)));
end

%% extract data
% GA
for n_subj_GA = 1:Num_GA
    tmp_idx = intersect(idx_subj_GA{n_subj_GA,1},idx_lateral_L);
    Fq_GA_L(n_subj_GA,:) = mean(Fq_All(tmp_idx,:));
    Pw_GA_L(n_subj_GA,:) = mean(Pw_All(tmp_idx,:));
    f_GA_L(n_subj_GA,:) = mean(f_All(tmp_idx,:));
    tmp_idx = intersect(idx_subj_GA{n_subj_GA,1},idx_lateral_R);
    Fq_GA_R(n_subj_GA,:) = mean(Fq_All(tmp_idx,:));
    Pw_GA_R(n_subj_GA,:) = mean(Pw_All(tmp_idx,:));
    f_GA_R(n_subj_GA,:) = mean(f_All(tmp_idx,:));
end
for n_subj_LA = 1:Num_LA
    tmp_idx = intersect(idx_subj_LA{n_subj_LA,1},idx_lateral_L);
    Fq_LA_L(n_subj_LA,:) = mean(Fq_All(tmp_idx,:));
    Pw_LA_L(n_subj_LA,:) = mean(Pw_All(tmp_idx,:));
    f_LA_L(n_subj_LA,:) = mean(f_All(tmp_idx,:));
    tmp_idx = intersect(idx_subj_LA{n_subj_LA,1},idx_lateral_R);
    Fq_LA_R(n_subj_LA,:) = mean(Fq_All(tmp_idx,:));
    Pw_LA_R(n_subj_LA,:) = mean(Pw_All(tmp_idx,:));
    f_LA_R(n_subj_LA,:) = mean(f_All(tmp_idx,:));
end
Fq_LA = (Fq_LA_L+Fq_LA_R)/2;
Fq_GA = (Fq_GA_L+Fq_GA_R)/2;
Pw_LA = (Pw_LA_L+Pw_LA_R)/2;
Pw_GA = (Pw_GA_L+Pw_GA_R)/2;
f_LA = (f_LA_L+f_LA_R)/2;
f_GA = (f_GA_L+f_GA_R)/2;

%% Delta(0.5-3Hz) TEST
delta_Fq_LA_L=mean(Fq_LA_L(:,6:31),2);delta_Fq_LA_R=mean(Fq_LA_R(:,6:31),2);
% [delta_Fq_LA_L, b_delta_Fq_LA_L, bint_delta_Fq_LA_L, ~] = regress_out(delta_Fq_LA_L1, age_LA);
% [delta_Fq_LA_R, b_delta_Fq_LA_R, bint_delta_Fq_LA_R, ~] = regress_out(delta_Fq_LA_R1, age_LA);
delta_Fq_GA_L=mean(Fq_GA_L(:,6:31),2);delta_Fq_GA_R=mean(Fq_GA_R(:,6:31),2);
% [delta_Fq_GA_L, b_delta_Fq_GA_L, bint_delta_Fq_GA_L, ~] = regress_out(delta_Fq_GA_L1, age_GA);
% [delta_Fq_GA_R, b_delta_Fq_GA_R, bint_delta_Fq_GA_R, ~] = regress_out(delta_Fq_GA_R1, age_GA);
delta_Fq_LA=mean(Fq_LA(:,6:31),2);delta_Fq_GA=mean(Fq_GA(:,6:31),2);
% [delta_Fq_LA, b_delta_Fq_LA, bint_delta_Fq_LA, ~] = regress_out(delta_Fq_LA1, age_LA);
% [delta_Fq_GA, b_delta_Fq_GA, bint_delta_Fq_GA, ~] = regress_out(delta_Fq_GA1, age_GA);
% delta_Pw_LA_L1=mean(Pw_LA_L(:,6:31),2);delta_Pw_LA_R1=mean(Pw_LA_R(:,6:31),2);
% [delta_Pw_LA_L, ~, ~] = regress_out(delta_Pw_LA_L1, old_LA);
% [delta_Pw_LA_R, ~, ~] = regress_out(delta_Pw_LA_R1, old_LA);
% delta_Pw_GA_L1=mean(Pw_GA_L(:,6:31),2);delta_Pw_GA_R1=mean(Pw_GA_R(:,6:31),2);
% [delta_Pw_GA_L, ~, ~] = regress_out(delta_Pw_GA_L1, old_GA);
% [delta_Pw_GA_R, ~, ~] = regress_out(delta_Pw_GA_R1, old_GA);
% delta_Pw_LA1=mean(Pw_LA(:,6:31),2);delta_Pw_GA1=mean(Pw_GA(:,6:31),2);
% [delta_Pw_LA, ~, ~] = regress_out(delta_Pw_LA1, old_LA);
% [delta_Pw_GA, ~, ~] = regress_out(delta_Pw_GA1, old_GA);

[h,p_d_L_v]=vartest2(delta_Fq_LA_L,delta_Fq_GA_L);
[h,p_d_R_v]=vartest2(delta_Fq_LA_R,delta_Fq_GA_R);
[h,p_d_v]=vartest2(delta_Fq_LA,delta_Fq_GA);

[tmp, p_delta_Fq_LA] = ttest(delta_Fq_LA_L,delta_Fq_LA_R);
[tmp, p_delta_Fq_GA] = ttest(delta_Fq_GA_L,delta_Fq_GA_R);
[tmp, p_delta_Fq_L] = ttest2(delta_Fq_LA_L,delta_Fq_GA_L);
[tmp, p_delta_Fq_R] = ttest2(delta_Fq_LA_R,delta_Fq_GA_R);
[tmp, p_delta_Fq] = ttest2(delta_Fq_LA,delta_Fq_GA);

% [tmp, p_delta_Pw_LA] = ttest(delta_Pw_LA_L,delta_Pw_LA_R);
% [tmp, p_delta_Pw_GA] = ttest(delta_Pw_GA_L,delta_Pw_GA_R);
% [tmp, p_delta_Pw_L] = ttest2(delta_Pw_LA_L,delta_Pw_GA_L);
% [tmp, p_delta_Pw_R] = ttest2(delta_Pw_LA_R,delta_Pw_GA_R);
% [tmp, p_delta_Pw] = ttest2(delta_Pw_LA,delta_Pw_GA); 

m_delta_Fq_LA = mean(delta_Fq_LA);sd_delta_Fq_LA = std(delta_Fq_LA);
m_delta_Fq_LA_L = mean(delta_Fq_LA_L);sd_delta_Fq_LA_L = std(delta_Fq_LA_L);
m_delta_Fq_LA_R = mean(delta_Fq_LA_R);sd_delta_Fq_LA_R = std(delta_Fq_LA_R);
m_delta_Fq_GA = mean(delta_Fq_GA);sd_delta_Fq_GA = std(delta_Fq_GA);
m_delta_Fq_GA_L = mean(delta_Fq_GA_L);sd_delta_Fq_GA_L = std(delta_Fq_GA_L);
m_delta_Fq_GA_R = mean(delta_Fq_GA_R);sd_delta_Fq_GA_R = std(delta_Fq_GA_R);
% m_delta_Pw_LA = mean(delta_Pw_LA);
% sd_delta_Pw_LA = std(delta_Pw_LA);
% m_delta_Pw_GA = mean(delta_Pw_GA);
% sd_delta_Pw_GA = std(delta_Pw_GA);

%% Theta(4-8Hz) TEST
theta_Fq_LA_L=mean(Fq_LA_L(:,41:81),2);theta_Fq_LA_R=mean(Fq_LA_R(:,41:81),2);
% [theta_Fq_LA_L, b_theta_Fq_LA_L, bint_theta_Fq_LA_L, ~] = regress_out(theta_Fq_LA_L1, age_LA);
% [theta_Fq_LA_R, b_theta_Fq_LA_R, bint_theta_Fq_LA_R, ~] = regress_out(theta_Fq_LA_R1, age_LA);
theta_Fq_GA_L=mean(Fq_GA_L(:,41:81),2);theta_Fq_GA_R=mean(Fq_GA_R(:,41:81),2);
% [theta_Fq_GA_L, b_theta_Fq_GA_L, bint_theta_Fq_GA_L, ~] = regress_out(theta_Fq_GA_L1, age_GA);
% [theta_Fq_GA_R, b_theta_Fq_GA_R, bint_theta_Fq_GA_R, ~] = regress_out(theta_Fq_GA_R1, age_GA);
theta_Fq_LA=mean(Fq_LA(:,41:81),2);theta_Fq_GA=mean(Fq_GA(:,41:81),2);
% [theta_Fq_LA, b_theta_Fq_LA, bint_theta_Fq_LA, ~] = regress_out(theta_Fq_LA1, age_LA);
% [theta_Fq_GA, b_theta_Fq_GA, bint_theta_Fq_GA, ~] = regress_out(theta_Fq_GA1, age_GA);
% theta_Pw_LA_L=mean(Pw_LA_L(:,41:81),2);theta_Pw_LA_R=mean(Pw_LA_R(:,41:81),2);
% theta_Pw_GA_L=mean(Pw_GA_L(:,41:81),2);theta_Pw_GA_R=mean(Pw_GA_R(:,41:81),2);
% theta_Pw_LA=mean(Pw_LA(:,41:81),2);theta_Pw_GA=mean(Pw_GA(:,41:81),2);

[h,p_t_L_v]=vartest2(theta_Fq_LA_L,theta_Fq_GA_L);
[h,p_t_R_v]=vartest2(theta_Fq_LA_R,theta_Fq_GA_R);
[h,p_t_v]=vartest2(theta_Fq_LA,theta_Fq_GA);

[tmp, p_theta_Fq_LA] = ttest(theta_Fq_LA_L,theta_Fq_LA_R);
[tmp, p_theta_Fq_GA] = ttest(theta_Fq_GA_L,theta_Fq_GA_R);
[tmp, p_theta_Fq_L] = ttest2(theta_Fq_LA_L,theta_Fq_GA_L);
[tmp, p_theta_Fq_R] = ttest2(theta_Fq_LA_R,theta_Fq_GA_R);
[tmp, p_theta_Fq] = ttest2(theta_Fq_LA,theta_Fq_GA);

% [tmp, p_theta_Pw_LA] = ttest(theta_Pw_LA_L,theta_Pw_LA_R);
% [tmp, p_theta_Pw_GA] = ttest(theta_Pw_GA_L,theta_Pw_GA_R);
% [tmp, p_theta_Pw_L] = ttest2(theta_Pw_LA_L,theta_Pw_GA_L);
% [tmp, p_theta_Pw_R] = ttest2(theta_Pw_LA_R,theta_Pw_GA_R);
% [tmp, p_theta_Pw] = ttest2(theta_Pw_LA,theta_Pw_GA); 

m_theta_Fq_LA = mean(theta_Fq_LA);sd_theta_Fq_LA = std(theta_Fq_LA);
m_theta_Fq_LA_L = mean(theta_Fq_LA_L);sd_theta_Fq_LA_L = std(theta_Fq_LA_L);
m_theta_Fq_LA_R = mean(theta_Fq_LA_R);sd_theta_Fq_LA_R = std(theta_Fq_LA_R);
m_theta_Fq_GA = mean(theta_Fq_GA);sd_theta_Fq_GA = std(theta_Fq_GA);
m_theta_Fq_GA_L = mean(theta_Fq_GA_L);sd_theta_Fq_GA_L = std(theta_Fq_GA_L);
m_theta_Fq_GA_R = mean(theta_Fq_GA_R);sd_theta_Fq_GA_R = std(theta_Fq_GA_R);
% m_theta_Pw_LA = mean(theta_Pw_LA);
% sd_theta_Pw_LA = std(theta_Pw_LA);
% m_theta_Pw_GA = mean(theta_Pw_GA);
% sd_theta_Pw_GA = std(theta_Pw_GA);

%% Alpha(8-14Hz) TEST
alpha_Fq_LA_L=mean(Fq_LA_L(:,81:141),2);alpha_Fq_LA_R=mean(Fq_LA_R(:,81:141),2);
% [alpha_Fq_LA_L, b_alpha_Fq_LA_L, bint_alpha_Fq_LA_LA, ~] = regress_out(alpha_Fq_LA_L1, age_LA);
% [alpha_Fq_LA_R, b_alpha_Fq_LA_R, bint_alpha_Fq_LA_R, ~] = regress_out(alpha_Fq_LA_R1, age_LA);
alpha_Fq_GA_L=mean(Fq_GA_L(:,81:141),2);alpha_Fq_GA_R=mean(Fq_GA_R(:,81:141),2);
% [alpha_Fq_GA_L, b_alpha_Fq_GA_L, bint_alpha_Fq_GA_L, ~] = regress_out(alpha_Fq_GA_L1, age_GA);
% [alpha_Fq_GA_R, b_alpha_Fq_GA_R, bint_alpha_Fq_GA_R, ~] = regress_out(alpha_Fq_GA_R1, age_GA);
alpha_Fq_LA=mean(Fq_LA(:,81:141),2);alpha_Fq_GA=mean(Fq_GA(:,81:141),2);
% [alpha_Fq_LA, b_alpha_Fq_LA, bint_alpha_Fq_LA, ~] = regress_out(alpha_Fq_LA1, age_LA);
% [alpha_Fq_GA, b_alpha_Fq_GA, bint_alpha_Fq_GA, ~] = regress_out(alpha_Fq_GA1, age_GA);
% alpha_Pw_LA_L=mean(Pw_LA_L(:,81:141),2);alpha_Pw_LA_R=mean(Pw_LA_R(:,81:141),2);
% alpha_Pw_GA_L=mean(Pw_GA_L(:,81:141),2);alpha_Pw_GA_R=mean(Pw_GA_R(:,81:141),2);
% alpha_Pw_LA=mean(Pw_LA(:,81:141),2);alpha_Pw_GA=mean(Pw_GA(:,81:141),2);

[h,p_a_L_v]=vartest2(alpha_Fq_LA_L,alpha_Fq_GA_L);
[h,p_a_R_v]=vartest2(alpha_Fq_LA_R,alpha_Fq_GA_R);
[h,p_a_v]=vartest2(alpha_Fq_LA,alpha_Fq_GA);

[tmp, p_alpha_Fq_LA] = ttest(alpha_Fq_LA_L,alpha_Fq_LA_R);
[tmp, p_alpha_Fq_GA] = ttest(alpha_Fq_GA_L,alpha_Fq_GA_R);
[tmp, p_alpha_Fq_L] = ttest2(alpha_Fq_LA_L,alpha_Fq_GA_L);
[tmp, p_alpha_Fq_R] = ttest2(alpha_Fq_LA_R,alpha_Fq_GA_R);
[tmp, p_alpha_Fq] = ttest2(alpha_Fq_LA,alpha_Fq_GA);

% [tmp, p_alpha_Pw_LA] = ttest(alpha_Pw_LA_L,alpha_Pw_LA_R);
% [tmp, p_alpha_Pw_GA] = ttest(alpha_Pw_GA_L,alpha_Pw_GA_R);
% [tmp, p_alpha_Pw_L] = ttest2(alpha_Pw_LA_L,alpha_Pw_GA_L);
% [tmp, p_alpha_Pw_R] = ttest2(alpha_Pw_LA_R,alpha_Pw_GA_R);
% [tmp, p_alpha_Pw] = ttest2(alpha_Pw_LA,alpha_Pw_GA); 

m_alpha_Fq_LA = mean(alpha_Fq_LA);sd_alpha_Fq_LA = std(alpha_Fq_LA);
m_alpha_Fq_LA_L = mean(alpha_Fq_LA_L);sd_alpha_Fq_LA_L = std(alpha_Fq_LA_L);
m_alpha_Fq_LA_R = mean(alpha_Fq_LA_R);sd_alpha_Fq_LA_R = std(alpha_Fq_LA_R);
m_alpha_Fq_GA = mean(alpha_Fq_GA);sd_alpha_Fq_GA = std(alpha_Fq_GA);
m_alpha_Fq_GA_L = mean(alpha_Fq_GA_L);sd_alpha_Fq_GA_L = std(alpha_Fq_GA_L);
m_alpha_Fq_GA_R = mean(alpha_Fq_GA_R);sd_alpha_Fq_GA_R = std(alpha_Fq_GA_R);
% m_alpha_Pw_LA = mean(alpha_Pw_LA);
% sd_alpha_Pw_LA = std(alpha_Pw_LA);
% m_alpha_Pw_GA = mean(alpha_Pw_GA);
% sd_alpha_Pw_GA = std(alpha_Pw_GA);

%% 0.5-25Hz All
all_Fq_LA_L=mean(Fq_LA_L(:,6:121),2);all_Fq_LA_R=mean(Fq_LA_R(:,6:121),2);
% [all_Fq_LA_L, b_all_Fq_LA_L, bint_all_Fq_LA_L, ~] = regress_out(all_Fq_LA_L1, age_LA);
% [all_Fq_LA_R, b_all_Fq_LA_R, bint_all_Fq_LA_R, ~] = regress_out(all_Fq_LA_R1, age_LA);
all_Fq_GA_L=mean(Fq_GA_L(:,6:121),2);all_Fq_GA_R=mean(Fq_GA_R(:,6:121),2);
% [all_Fq_GA_L, b_all_Fq_GA_L, bint_all_Fq_GA_L, ~] = regress_out(all_Fq_GA_L1, age_GA);
% [all_Fq_GA_R, b_all_Fq_GA_R, bint_all_Fq_GA_R, ~] = regress_out(all_Fq_GA_R1, age_GA);
all_Fq_LA=mean(Fq_LA(:,6:121),2);all_Fq_GA=mean(Fq_GA(:,6:121),2);
% [all_Fq_LA, b_all_Fq_LA, bint_all_Fq_LA, ~] = regress_out(all_Fq_LA1, age_LA);
% [all_Fq_GA, b_all_Fq_GA, bint_all_Fq_GA, ~] = regress_out(all_Fq_GA1, age_GA);

[h,p_a_L_v]=vartest2(alpha_Fq_LA_L,alpha_Fq_GA_L);
[h,p_a_R_v]=vartest2(alpha_Fq_LA_R,alpha_Fq_GA_R);
[h,p_a_v]=vartest2(alpha_Fq_LA,alpha_Fq_GA);

[tmp, p_all_Fq_LA] = ttest(all_Fq_LA_L,all_Fq_LA_R);
[tmp, p_all_Fq_GA] = ttest(all_Fq_GA_L,all_Fq_GA_R);
[tmp, p_all_Fq_L] = ttest2(all_Fq_LA_L,all_Fq_GA_L);
[tmp, p_all_Fq_R] = ttest2(all_Fq_LA_R,all_Fq_GA_R);
[tmp, p_all_Fq] = ttest2(all_Fq_LA,all_Fq_GA);

m_all_Fq_LA = mean(all_Fq_LA);sd_all_Fq_LA = std(all_Fq_LA);
m_all_Fq_LA_L = mean(all_Fq_LA_L);sd_all_Fq_LA_L = std(all_Fq_LA_L);
m_all_Fq_LA_R = mean(all_Fq_LA_R);sd_all_Fq_LA_R = std(all_Fq_LA_R);
m_all_Fq_GA = mean(all_Fq_GA);sd_all_Fq_GA = std(all_Fq_GA);
m_all_Fq_GA_L = mean(all_Fq_GA_L);sd_all_Fq_GA_L = std(all_Fq_GA_L);
m_all_Fq_GA_R = mean(all_Fq_GA_R);sd_all_Fq_GA_R = std(all_Fq_GA_R);

%% Beta(14-31Hz) TEST
beta_Fq_LA_L=mean(Fq_LA_L(:,141:311),2);beta_Fq_LA_R=mean(Fq_LA_R(:,141:311),2);
% [beta_Fq_LA_L, b_beta_Fq_LA_L, bint_beta_Fq_LA_L, ~] = regress_out(beta_Fq_LA_L1, age_LA);
% [beta_Fq_LA_R, b_beta_Fq_LA_R, bint_beta_Fq_LA_R, ~] = regress_out(beta_Fq_LA_R1, age_LA);
beta_Fq_GA_L=mean(Fq_GA_L(:,141:311),2);beta_Fq_GA_R=mean(Fq_GA_R(:,141:311),2);
% [beta_Fq_GA_L, b_beta_Fq_GA_L, bint_beta_Fq_GA_L, ~] = regress_out(beta_Fq_GA_L1, age_GA);
% [beta_Fq_GA_R, b_beta_Fq_GA_R, bint_beta_Fq_GA_R, ~] = regress_out(beta_Fq_GA_R1, age_GA);
beta_Fq_LA=mean(Fq_LA(:,141:311),2);beta_Fq_GA=mean(Fq_GA(:,141:311),2);
% [beta_Fq_LA, b_beta_Fq_LA, bint_beta_Fq_LA, ~] = regress_out(beta_Fq_LA1, age_LA);
% [beta_Fq_GA, b_beta_Fq_GA, bint_beta_Fq_GA, ~] = regress_out(beta_Fq_GA1, age_GA);
% beta_Pw_LA_L=mean(Pw_LA_L(:,141:311),2);beta_Pw_LA_R=mean(Pw_LA_R(:,141:311),2);
% beta_Pw_GA_L=mean(Pw_GA_L(:,141:311),2);beta_Pw_GA_R=mean(Pw_GA_R(:,141:311),2);
% beta_Pw_LA=mean(Pw_LA(:,141:311),2);beta_Pw_GA=mean(Pw_GA(:,141:311),2);

[h,p_b_L_v]=vartest2(beta_Fq_LA_L,beta_Fq_GA_L);
[h,p_b_R_v]=vartest2(beta_Fq_LA_R,beta_Fq_GA_R);
[h,p_b_v]=vartest2(beta_Fq_LA,beta_Fq_GA);

[tmp, p_beta_Fq_LA] = ttest(beta_Fq_LA_L,beta_Fq_LA_R);
[tmp, p_beta_Fq_GA] = ttest(beta_Fq_GA_L,beta_Fq_GA_R);
[tmp, p_beta_Fq_L] = ttest2(beta_Fq_LA_L,beta_Fq_GA_L);
[tmp, p_beta_Fq_R] = ttest2(beta_Fq_LA_R,beta_Fq_GA_R);
[tmp, p_beta_Fq] = ttest2(beta_Fq_LA,beta_Fq_GA);

% [tmp, p_beta_Pw_LA] = ttest(beta_Pw_LA_L,beta_Pw_LA_R);
% [tmp, p_beta_Pw_GA] = ttest(beta_Pw_GA_L,beta_Pw_GA_R);
% [tmp, p_beta_Pw_L] = ttest2(beta_Pw_LA_L,beta_Pw_GA_L);
% [tmp, p_beta_Pw_R] = ttest2(beta_Pw_LA_R,beta_Pw_GA_R);
% [tmp, p_beta_Pw] = ttest2(beta_Pw_LA,beta_Pw_GA); 

m_beta_Fq_LA = mean(beta_Fq_LA);sd_beta_Fq_LA = std(beta_Fq_LA);
m_beta_Fq_LA_L = mean(beta_Fq_LA_L);sd_beta_Fq_LA_L = std(beta_Fq_LA_L);
m_beta_Fq_LA_R = mean(beta_Fq_LA_R);sd_beta_Fq_LA_R = std(beta_Fq_LA_R);
m_beta_Fq_GA = mean(beta_Fq_GA);sd_beta_Fq_GA = std(beta_Fq_GA);
m_beta_Fq_GA_L = mean(beta_Fq_GA_L);sd_beta_Fq_GA_L = std(beta_Fq_GA_L);
m_beta_Fq_GA_R = mean(beta_Fq_GA_R);sd_beta_Fq_GA_R = std(beta_Fq_GA_R);
% m_beta_Pw_LA = mean(beta_Pw_LA);
% sd_beta_Pw_LA = std(beta_Pw_LA);
% m_beta_Pw_GA = mean(beta_Pw_GA);
% sd_beta_Pw_GA = std(beta_Pw_GA);

%% Gamma(31-100Hz) TEST
gamma_Fq_LA_L=mean(Fq_LA_L(:,311:1001),2);gamma_Fq_LA_R=mean(Fq_LA_R(:,311:1001),2);
% [gamma_Fq_LA_L, b_gamma_Fq_LA_L, bint_gamma_Fq_LA_L, ~] = regress_out(gamma_Fq_LA_L1, age_LA);
% [gamma_Fq_LA_R, b_gamma_Fq_LA_R, bint_gamma_Fq_LA_R, ~] = regress_out(gamma_Fq_LA_R1, age_LA);
gamma_Fq_GA_L=mean(Fq_GA_L(:,311:1001),2);gamma_Fq_GA_R=mean(Fq_GA_R(:,311:1001),2);
% [gamma_Fq_GA_L, b_gamma_Fq_GA_L, bint_gamma_Fq_GA_L, ~] = regress_out(gamma_Fq_GA_L1, age_GA);
% [gamma_Fq_GA_R, b_gamma_Fq_GA_R, bint_gamma_Fq_GA_R, ~] = regress_out(gamma_Fq_GA_R1, age_GA);
gamma_Fq_LA=mean(Fq_LA(:,311:1001),2);gamma_Fq_GA=mean(Fq_GA(:,311:1001),2);
% [gamma_Fq_LA, b_gamma_Fq_LA, bint_gamma_Fq_LA, ~] = regress_out(gamma_Fq_LA1, age_LA);
% [gamma_Fq_GA, b_gamma_Fq_GA, bint_gamma_Fq_GA, ~] = regress_out(gamma_Fq_GA1, age_GA);
% gamma_Pw_LA_L=mean(Pw_LA_L(:,311:1001),2);gamma_Pw_LA_R=mean(Pw_LA_R(:,311:1001),2);
% gamma_Pw_GA_L=mean(Pw_GA_L(:,311:1001),2);gamma_Pw_GA_R=mean(Pw_GA_R(:,311:1001),2);
% gamma_Pw_LA=mean(Pw_LA(:,311:1001),2);gamma_Pw_GA=mean(Pw_GA(:,311:1001),2);

[h,p_g_L_v]=vartest2(gamma_Fq_LA_L,gamma_Fq_GA_L);
[h,p_g_R_v]=vartest2(gamma_Fq_LA_R,gamma_Fq_GA_R);
[h,p_g_v]=vartest2(gamma_Fq_LA,gamma_Fq_GA);

[tmp, p_gamma_Fq_LA] = ttest(gamma_Fq_LA_L,gamma_Fq_LA_R);
[tmp, p_gamma_Fq_GA] = ttest(gamma_Fq_GA_L,gamma_Fq_GA_R);
[tmp, p_gamma_Fq_L] = ttest2(gamma_Fq_LA_L,gamma_Fq_GA_L);
[tmp, p_gamma_Fq_R] = ttest2(gamma_Fq_LA_R,gamma_Fq_GA_R);
[tmp, p_gamma_Fq] = ttest2(gamma_Fq_LA,gamma_Fq_GA);

% [tmp, p_gamma_Pw_LA] = ttest(gamma_Pw_LA_L,gamma_Pw_LA_R);
% [tmp, p_gamma_Pw_GA] = ttest(gamma_Pw_GA_L,gamma_Pw_GA_R);
% [tmp, p_gamma_Pw_L] = ttest2(gamma_Pw_LA_L,gamma_Pw_GA_L);
% [tmp, p_gamma_Pw_R] = ttest2(gamma_Pw_LA_R,gamma_Pw_GA_R);
% [tmp, p_gamma_Pw] = ttest2(gamma_Pw_LA,gamma_Pw_GA); 

m_gamma_Fq_LA = mean(gamma_Fq_LA);sd_gamma_Fq_LA = std(gamma_Fq_LA);
m_gamma_Fq_LA_L = mean(gamma_Fq_LA_L);sd_gamma_Fq_LA_L = std(gamma_Fq_LA_L);
m_gamma_Fq_LA_R = mean(gamma_Fq_LA_R);sd_gamma_Fq_LA_R = std(gamma_Fq_LA_R);
m_gamma_Fq_GA = mean(gamma_Fq_GA);sd_gamma_Fq_GA = std(gamma_Fq_GA);
m_gamma_Fq_GA_L = mean(gamma_Fq_GA_L);sd_gamma_Fq_GA_L = std(gamma_Fq_GA_L);
m_gamma_Fq_GA_R = mean(gamma_Fq_GA_R);sd_gamma_Fq_GA_R = std(gamma_Fq_GA_R);

% m_gamma_Pw_LA = mean(gamma_Pw_LA
% sd_gamma_Pw_LA = std(gamma_Pw_LA);
% m_gamma_Pw_GA = mean(gamma_Pw_GA);
% sd_gamma_Pw_GA = std(gamma_Pw_GA);

%% 2-way anova
%delta
delta_L=[delta_Fq_LA_L;delta_Fq_GA_L];
delta_R=[delta_Fq_LA_R;delta_Fq_GA_R];
anesthesia=[zeros(n_subj_LA,1)+1;zeros(n_subj_GA,1)+2];
d_t=table(anesthesia,delta_L,delta_R,'VariableNames',{'anesthesia','L','R'});
d_meas=table([1 2]','VariableNames',{'STN_side'});
d_rm=fitrm(d_t,'L,R ~ anesthesia','WithinDesign',d_meas);
d_ranovatbl=ranova(d_rm)
d_anovatbl=anova(d_rm)
[tmp, p_delta_LR] = ttest2(delta_L,delta_R);

%theta
theta_L=[theta_Fq_LA_L;theta_Fq_GA_L];
theta_R=[theta_Fq_LA_R;theta_Fq_GA_R];
anesthesia=[zeros(n_subj_LA,1)+1;zeros(n_subj_GA,1)+2];
t_t=table(anesthesia,theta_L,theta_R,'VariableNames',{'anesthesia','L','R'});
t_meas=table([1 2]','VariableNames',{'STN_side'});
t_rm=fitrm(t_t,'L,R ~ anesthesia','WithinDesign',t_meas);
t_ranovatbl=ranova(t_rm)
t_anovatbl=anova(t_rm)
[tmp, p_theta_LR] = ttest2(theta_L,theta_R);

%alpha
alpha_L=[alpha_Fq_LA_L;alpha_Fq_GA_L];
alpha_R=[alpha_Fq_LA_R;alpha_Fq_GA_R];
anesthesia=[zeros(n_subj_LA,1)+1;zeros(n_subj_GA,1)+2];
a_t=table(anesthesia,alpha_L,alpha_R,'VariableNames',{'anesthesia','L','R'});
a_meas=table([1 2]','VariableNames',{'STN_side'});
a_rm=fitrm(a_t,'L,R ~ anesthesia','WithinDesign',a_meas);
a_ranovatbl=ranova(a_rm)
a_anovatbl=anova(a_rm)
[tmp, p_alpha_LR] = ttest2(alpha_L,alpha_R);

%all

all_L=[all_Fq_LA_L;all_Fq_GA_L];
all_R=[all_Fq_LA_R;all_Fq_GA_R];
anesthesia=[zeros(n_subj_LA,1)+1;zeros(n_subj_GA,1)+2];
all_t=table(anesthesia,all_L,all_R,'VariableNames',{'anesthesia','L','R'});
all_meas=table([1 2]','VariableNames',{'STN_side'});
all_rm=fitrm(all_t,'L,R ~ anesthesia','WithinDesign',all_meas);
all_ranovatbl=ranova(all_rm);
all_anovatbl=anova(all_rm)


%beta
beta_L=[beta_Fq_LA_L;beta_Fq_GA_L];
beta_R=[beta_Fq_LA_R;beta_Fq_GA_R];
anesthesia=[zeros(n_subj_LA,1)+1;zeros(n_subj_GA,1)+2];
b_t=table(anesthesia,beta_L,beta_R,'VariableNames',{'anesthesia','L','R'});
b_meas=table([1 2]','VariableNames',{'STN_side'});
b_rm=fitrm(b_t,'L,R ~ anesthesia','WithinDesign',b_meas);
b_ranovatbl=ranova(b_rm)
b_anovatbl=anova(b_rm)
[tmp, p_beta_LR] = ttest2(beta_L,beta_R);

%gamma
gamma_L=[gamma_Fq_LA_L;gamma_Fq_GA_L];
gamma_R=[gamma_Fq_LA_R;gamma_Fq_GA_R];
anesthesia=[zeros(n_subj_LA,1)+1;zeros(n_subj_GA,1)+2];
g_t=table(anesthesia,gamma_L,gamma_R,'VariableNames',{'anesthesia','L','R'});
g_meas=table([1 2]','VariableNames',{'STN_side'});
g_rm=fitrm(g_t,'L,R ~ anesthesia','WithinDesign',g_meas);
g_ranovatbl=ranova(g_rm)
g_anovatbl=anova(g_rm)
[tmp, p_gamma_LR] = ttest2(gamma_L,gamma_R);