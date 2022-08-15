clear all

%% Load Results
dirname = 'Results\';
w = what(dirname);
resultfilenames = w.mat;

%% extract index for matrix concatenation
Num_LA = 0; Num_GA = 0;
str_subject_LA = {}; str_subject_GA = {};

for n_file=1:numel(resultfilenames) % START
    n_file
    filename = [dirname,resultfilenames{n_file}];
    load (filename)
    MBI_All(n_file,1) = MBI;
    FR_All(n_file,1) = FR;

    IDX_Group_All{n_file,1} = resultfilenames{n_file}(1:2);
    IDX_Subject_All{n_file,1} = resultfilenames{n_file}(4:5);
    IDX_Lateral_All{n_file,1} = resultfilenames{n_file}(7);
    IDX_Depth_All{n_file,1} = resultfilenames{n_file}(9:12);
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
    MBI_GA_L(n_subj_GA,1) = nanmean(MBI_All(tmp_idx));
    FR_GA_L(n_subj_GA,1) = nanmean(FR_All(tmp_idx));
    tmp_idx = intersect(idx_subj_GA{n_subj_GA,1},idx_lateral_R);
    MBI_GA_R(n_subj_GA,1) = nanmean(MBI_All(tmp_idx));
    FR_GA_R(n_subj_GA,1) = nanmean(FR_All(tmp_idx));
end
for n_subj_LA = 1:Num_LA
    tmp_idx = intersect(idx_subj_LA{n_subj_LA,1},idx_lateral_L);
    MBI_LA_L(n_subj_LA,1) = nanmean(MBI_All(tmp_idx));
    FR_LA_L(n_subj_LA,1) = nanmean(FR_All(tmp_idx));
    tmp_idx = intersect(idx_subj_LA{n_subj_LA,1},idx_lateral_R);
    MBI_LA_R(n_subj_LA,1) = nanmean(MBI_All(tmp_idx));
    FR_LA_R(n_subj_LA,1) = nanmean(FR_All(tmp_idx));
end
MBI_LA = (MBI_LA_L+MBI_LA_R)/2;
MBI_GA = (MBI_GA_L+MBI_GA_R)/2;
FR_LA = (FR_LA_L+FR_LA_R)/2;
FR_GA = (FR_GA_L+FR_GA_R)/2;
%% TEST
[tmp, p_MBI_LA] = ttest(MBI_LA_L,MBI_LA_R);
[tmp, p_MBI_GA] = ttest(MBI_GA_L,MBI_GA_R);
[tmp, p_MBI_L] = ttest2(MBI_LA_L,MBI_GA_L);
[tmp, p_MBI_R] = ttest2(MBI_LA_R,MBI_GA_R);
[tmp, p_MBI] = ttest2(MBI_LA,MBI_GA);

[tmp, p_FR_LA] = ttest(FR_LA_L,FR_LA_R,0.05);
[tmp, p_FR_GA] = ttest(FR_GA_L,FR_GA_R,0.05);
[tmp, p_FR_L] = ttest2(FR_LA_L,FR_GA_L,0.05);
[tmp, p_FR_R] = ttest2(FR_LA_R,FR_GA_R,0.05);
[tmp, p_FR] = ttest2(FR_LA,FR_GA,0.05); 
% tmp=0,零假设在5%的置信度下不被拒绝,两者在统计上可看作统一分布数据;
% tmp=1,零假设被拒绝,两者来自不同分布

m_MBI_LA = mean(MBI_LA);sd_MBI_LA = std(MBI_LA);
m_MBI_LA_L = mean(MBI_LA_L);sd_MBI_LA_L = std(MBI_LA_L);
m_MBI_LA_R = mean(MBI_LA_R);sd_MBI_LA_R = std(MBI_LA_R);
m_MBI_GA = mean(MBI_GA);sd_MBI_GA = std(MBI_GA);
m_MBI_GA_L = mean(MBI_GA_L);sd_MBI_GA_L = std(MBI_GA_L);
m_MBI_GA_R = mean(MBI_GA_R);sd_MBI_GA_R = std(MBI_GA_R);
m_FR_LA = mean(FR_LA);sd_FR_LA = std(FR_LA);
m_FR_LA_L = mean(FR_LA_L);sd_FR_LA_L = std(FR_LA_L);
m_FR_LA_R = mean(FR_LA_R);sd_FR_LA_R = std(FR_LA_R);
m_FR_GA = mean(FR_GA);sd_FR_GA = std(FR_GA);
m_FR_GA_L = mean(FR_GA_L);sd_FR_GA_L = std(FR_GA_L);
m_FR_GA_R = mean(FR_GA_R);sd_FR_GA_R = std(FR_GA_R);

%% 2-way anova
%MBI
MBI_L=[MBI_LA_L;MBI_GA_L];
MBI_R=[MBI_LA_R;MBI_GA_R];
anesthesia=[zeros(n_subj_LA,1)+1;zeros(n_subj_GA,1)+2];
MBI_t=table(anesthesia,MBI_L,MBI_R,'VariableNames',{'anesthesia','L','R'});
MBI_meas=table([1 2]','VariableNames',{'STN_side'});
MBI_rm=fitrm(MBI_t,'L,R ~ anesthesia','WithinDesign',MBI_meas);
MBI_ranovatbl=ranova(MBI_rm)
MBI_anovatbl=anova(MBI_rm)
[tmp, p_MBI_LR] = ttest2(MBI_L,MBI_R); 

%FR
FR_L=[FR_LA_L;FR_GA_L];
FR_R=[FR_LA_R;FR_GA_R];
anesthesia=[zeros(n_subj_LA,1)+1;zeros(n_subj_GA,1)+2];
FR_t=table(anesthesia,FR_L,FR_R,'VariableNames',{'anesthesia','L','R'});
FR_meas=table([1 2]','VariableNames',{'STN_side'});
FR_rm=fitrm(FR_t,'L,R ~ anesthesia','WithinDesign',FR_meas);
FR_ranovatbl=ranova(FR_rm)
FR_anovatbl=anova(FR_rm)
[tmp, p_FR_LR] = ttest2(FR_L,FR_R); 