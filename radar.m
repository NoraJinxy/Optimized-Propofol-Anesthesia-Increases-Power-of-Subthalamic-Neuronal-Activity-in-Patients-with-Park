clear all
close all
clc
%{
awake_baseline=[64.5,21,837.5,19.08,95.83,7.08,31.32,41.67,6.5,13.75,27.91,22.82];
awake_6months=[44,16.5,588,15.58,114,5.58,20.28,18.75,5.5,7.67,27.46,23.73];
TCI_baseline=[55,21,800,9.4,108.9,9.7,32.79,29.3,4.5,6.7,28.5,24.2];
TCI_6months=[30,19,375,6.7,116.85,8.7,15.39,17.1,2,6.5,27.9,25.3];
range=[0,1.2];
labels={'MDS-UPDRS III off','MDS-UPDRS III on','LEDD','PSQI','PDSS','ESS','PDQ-39','NMSS','HAMD','HAMA','MMSE','MOCA'};
draw_radar1(awake_baseline,awake_6months,TCI_baseline,TCI_6months,range,labels)
%}
%% Ô­
%{
awake_change=[0.33,-0.06,0.24,0.06,0.14,0.15,0.34,0.57,0.46,0.42,-0.07,0.05];
TCI_cahnge=[0.55,0.21,0.43,0.36,0.08,0.01,0.46,0.49,0.17,0.18,0,0.07];
range=[0,1.6];
labels={'MDS-UPDRS III off','MDS-UPDRS III on','LEDD','PSQI','PDSS','ESS','PDQ-39','NMSS','HAMD','HAMA','MMSE','MOCA'};
draw_radar1(awake_change,TCI_cahnge,range,labels)
%}
%% 2021/02/21¸Ä£¬È¥µôLEDD
awake_change=[0.33,-0.06,0.06,0.14,0.15,0.34,0.57,0.46,0.42,-0.07,0.05];
TCI_cahnge=[0.55,0.21,0.36,0.08,0.01,0.46,0.49,0.17,0.18,0,0.07];
range=[0,1.6];
labels={'MDS-UPDRS III off','MDS-UPDRS III on','PSQI','PDSS','ESS','PDQ-39','NMSS','HAMD','HAMA','MMSE','MOCA'};
draw_radar1(awake_change,TCI_cahnge,range,labels)