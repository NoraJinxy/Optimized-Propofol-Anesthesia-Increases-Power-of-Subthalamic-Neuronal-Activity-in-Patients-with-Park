clc;clear all;close all
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

delet_LA=find(score_pre_LA==0);
delet_GA=find(score_pre_GA==0);

[tmp, p_age] = ttest2(age_LA,age_GA);
[tmp, p_score_pre] = ttest2(score_pre_GA,score_pre_LA);
[tmp, p_score_post] = ttest2(score_post_GA,score_post_LA);