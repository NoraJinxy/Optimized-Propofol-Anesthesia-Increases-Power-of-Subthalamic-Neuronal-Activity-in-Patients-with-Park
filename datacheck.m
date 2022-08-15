clear all

load datainfo.mat % extracted from the EXCEL file
load PARAM.mat % extracted from text files
LA_L_dir_subj = LA(:,2);
LA_R_dir_subj = LA(:,5);
GA_L_dir_subj = GA(:,2);
GA_R_dir_subj = GA(:,5);

%% extract new parameters from each file & rename files
N_Files = numel(matfiles);
for n_file=1:N_Files
    % group
    if strcmp(PARAM1_ALL{n_file}(1),'A')
        NewParam{n_file,1}.group = 'LA';
        NewParamArray.group{n_file,1} = 'LA';
    elseif strcmp(PARAM1_ALL{n_file}(1),'B')
        NewParam{n_file,1}.group = 'GA';
        NewParamArray.group{n_file,1} = 'GA';
    end
    
    % lateral
	if strcmp(PARAM3_ALL{n_file},'17')
        NewParam{n_file}.lateral = 'L';
        NewParamArray.lateral{n_file,1} = 'L';
    elseif strcmp(PARAM3_ALL{n_file},'18')
        NewParam{n_file}.lateral = 'R';
        NewParamArray.lateral{n_file,1} = 'R';
    end   
    % subject
    dir_i = str2double(PARAM2_ALL{n_file});
    if strcmp(NewParam{n_file,1}.group,'LA')&&strcmp(PARAM3_ALL{n_file},'17')
        NewParam{n_file}.subject = find(LA_L_dir_subj==dir_i);
        NewParamArray.subject(n_file,1) = find(LA_L_dir_subj==dir_i);
    elseif strcmp(NewParam{n_file,1}.group,'LA')&&strcmp(PARAM3_ALL{n_file},'18')
        NewParam{n_file}.subject = find(LA_R_dir_subj==dir_i);
        NewParamArray.subject(n_file,1) = find(LA_R_dir_subj==dir_i);
    elseif strcmp(NewParam{n_file,1}.group,'GA')&&strcmp(PARAM3_ALL{n_file},'17')
        NewParam{n_file}.subject = find(GA_L_dir_subj==dir_i);
        NewParamArray.subject(n_file,1) = find(GA_L_dir_subj==dir_i);
    elseif strcmp(NewParam{n_file,1}.group,'GA')&&strcmp(PARAM3_ALL{n_file},'18')
        NewParam{n_file}.subject = find(GA_R_dir_subj==dir_i);
        NewParamArray.subject(n_file,1) = find(GA_R_dir_subj==dir_i);
    end
	% depth
    NewParam{n_file}.depth = str2double(PARAM7_ALL{n_file}(1:4))/10;
    NewParamArray.depth(n_file,1) = str2double(PARAM7_ALL{n_file}(1:4))/10;
    NewParam{n_file}.depth_str = PARAM7_ALL{n_file};
    NewParamArray.depth_str{n_file,1} = PARAM7_ALL{n_file};
    % other strings
    NewParam{n_file}.others = PARAM6_ALL{n_file};
    NewParamArray.others{n_file,1} = PARAM6_ALL{n_file};
    % New File Name
    newfilename{n_file,1} = [NewParam{n_file}.group,'_',num2str(NewParam{n_file}.subject,'%02.0f\n'),'_',NewParam{n_file}.lateral,'_',NewParam{n_file}.depth_str,'_',NewParam{n_file}.others];
end

%% extract info for each subject
idx_group_LA = find(strcmp(NewParamArray.group,'LA'));
idx_group_GA = find(strcmp(NewParamArray.group,'GA'));
idx_lateral_L = find(strcmp(NewParamArray.lateral,'L'));
idx_lateral_R = find(strcmp(NewParamArray.lateral,'R'));
idx_LA_L = intersect(idx_group_LA,idx_lateral_L);
idx_LA_R = intersect(idx_group_LA,idx_lateral_R);
idx_GA_L = intersect(idx_group_GA,idx_lateral_L);
idx_GA_R = intersect(idx_group_GA,idx_lateral_R);
for n_subject=1:max(NewParamArray.subject)
    % index
    idx_subject = find(NewParamArray.subject==n_subject);
    IDX_files{1,n_subject,1} = intersect(idx_LA_L,idx_subject);
    IDX_files{1,n_subject,2} = intersect(idx_LA_R,idx_subject);
    IDX_files{2,n_subject,1} = intersect(idx_GA_L,idx_subject);
    IDX_files{2,n_subject,2} = intersect(idx_GA_R,idx_subject);
    % number of files
    for n_group=1:2
        for n_lateral=1:2
            NUM_files(n_group,n_subject,n_lateral) = numel(IDX_files{n_group,n_subject,n_lateral});
            idx_tmp = IDX_files{n_group,n_subject,n_lateral};
            depth_range = NewParamArray.depth(idx_tmp,1);
            NUM_Depth(n_group,n_subject,n_lateral)=numel(unique(depth_range));
            MIN_Depth(n_group,n_subject,n_lateral)=min(depth_range);
            MAX_Depth(n_group,n_subject,n_lateral)=max(depth_range);
        end
    end
end
NUM_files_LA_L = squeeze(NUM_files(1,:,1))';
NUM_files_LA_R = squeeze(NUM_files(1,:,2))';
NUM_files_GA_L = squeeze(NUM_files(2,:,1))';
NUM_files_GA_R = squeeze(NUM_files(2,:,2))';
%sum(NUM_files_LA_L)+sum(NUM_files_LA_R)+sum(NUM_files_GA_L)+sum(NUM_files_GA_R)
NUM_Depth_LA_L = squeeze(NUM_Depth(1,:,1))';
NUM_Depth_LA_R = squeeze(NUM_Depth(1,:,2))';
NUM_Depth_GA_L = squeeze(NUM_Depth(2,:,1))';
NUM_Depth_GA_R = squeeze(NUM_Depth(2,:,2))';
MIN_Depth_LA_L = squeeze(MIN_Depth(1,:,1))';
MIN_Depth_LA_R = squeeze(MIN_Depth(1,:,2))';
MIN_Depth_GA_L = squeeze(MIN_Depth(2,:,1))';
MIN_Depth_GA_R = squeeze(MIN_Depth(2,:,2))';
MAX_Depth_LA_L = squeeze(MAX_Depth(1,:,1))';
MAX_Depth_LA_R = squeeze(MAX_Depth(1,:,2))';
MAX_Depth_GA_L = squeeze(MAX_Depth(2,:,1))';
MAX_Depth_GA_R = squeeze(MAX_Depth(2,:,2))';

%% show depth
close all
figure('unit','normalized','position',[0    0.06    1.0000    0.83])
subplot(2,1,1)
hold on; box on;
for n_subject=1:60
    scatter((n_subject-0.1)*ones(NUM_files_LA_L(n_subject),1),NewParamArray.depth(IDX_files{1,n_subject,1}),'b.');
    scatter((n_subject+0.1)*ones(NUM_files_LA_R(n_subject),1),NewParamArray.depth(IDX_files{1,n_subject,2}),'r.');
end
set(gca,'ygrid','on','xlim',[0 61],'ylim',[-8 8])
legend('left','right')
xlabel('# subject'); ylabel('Depth');
title('LA Group','fontsize',12)
subplot(2,1,2)
hold on; box on;
for n_subject=1:60
    scatter((n_subject-0.1)*ones(NUM_files_GA_L(n_subject),1),NewParamArray.depth(IDX_files{2,n_subject,1}),'b.');
    scatter((n_subject+0.1)*ones(NUM_files_GA_R(n_subject),1),NewParamArray.depth(IDX_files{2,n_subject,2}),'r.');
end
set(gca,'ygrid','on','xlim',[0 61],'ylim',[-8 8])
legend('left','right')
xlabel('# subject'); ylabel('Depth');
title('GA Group','fontsize',12)

save ('NewParam.mat','newfilename','NUM_*','MAX_*','MIN_*','NewParam','NewParamArray')

%% rename files
for n_file=1:N_Files
    n_file
	loadfile = ['MATFiles\',matfiles{n_file}];
    load (loadfile)
    param = NewParam{n_file};
    savefile = ['MATFilesNew\',newfilename{n_file}];
    save (savefile,'x','param');
end