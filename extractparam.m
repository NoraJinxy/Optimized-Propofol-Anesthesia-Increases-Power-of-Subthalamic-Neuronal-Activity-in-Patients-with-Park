clear all

s = what('MATFiles');
matfiles = s.mat;
N_Files = numel(matfiles);

for n_file=1:N_Files
    count_ = 0;
    for i=1:numel(matfiles{n_file})
        strinfo = matfiles{n_file};
        if strcmp(strinfo(i),'_')
            count_ = count_+1;
            idx_(count_) = i;
        end
    end
    
    PARAM1_ALL{n_file,1} = strinfo(1:idx_(1)-1); % A/B?
    PARAM2_ALL{n_file,1} = strinfo(idx_(1)+1:idx_(2)-1); % subject/dir
    PARAM3_ALL{n_file,1} = strinfo(idx_(2)+1:idx_(3)-1); % left or right
    PARAM4_ALL{n_file,1} = strinfo(idx_(3)+1:idx_(4)-1); % ?
    PARAM5_ALL{n_file,1} = strinfo(idx_(4)+1:idx_(5)-1); % ?
    PARAM6_ALL{n_file,1} = strinfo(idx_(5)+1:idx_(6)-1); % ?
    PARAM7_ALL{n_file,1} = strinfo(idx_(6)+1:end-4); % depth
end
PARAM1_UNQ = unique(PARAM1_ALL);
PARAM2_UNQ = unique(PARAM2_ALL);
PARAM3_UNQ = unique(PARAM3_ALL);
PARAM4_UNQ = unique(PARAM4_ALL);
PARAM5_UNQ = unique(PARAM5_ALL);
PARAM6_UNQ = unique(PARAM6_ALL);
PARAM7_UNQ = unique(PARAM7_ALL);

save ('PARAM.mat','PARAM*','matfiles')