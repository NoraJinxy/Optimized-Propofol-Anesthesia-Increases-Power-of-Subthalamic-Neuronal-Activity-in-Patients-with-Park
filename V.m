close all
clc
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