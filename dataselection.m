% select files
clear all

load('NewParam.mat')

N_Selected = 8; % number of selected files for each subject at each side

%%
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
            Range_Depth{n_group,n_subject,n_lateral} = unique(depth_range);
            NUM_Depth(n_group,n_subject,n_lateral)=numel(unique(depth_range));
            for nn=1:NUM_Depth(n_group,n_subject,n_lateral)
                Idx_EachDepth{n_group,n_subject,n_lateral}{nn,1} = find(depth_range==Range_Depth{n_group,n_subject,n_lateral}(nn));
                Num_EachDepth{n_group,n_subject,n_lateral}(nn,1) = numel(Idx_EachDepth{n_group,n_subject,n_lateral}{nn});
            end
            MIN_Depth(n_group,n_subject,n_lateral)=min(depth_range);
            MAX_Depth(n_group,n_subject,n_lateral)=max(depth_range);
            % remove some files
            n_remove_idx(n_group,n_subject,n_lateral) = NUM_Depth(n_group,n_subject,n_lateral) - N_Selected;
            if rem(n_remove_idx(n_group,n_subject,n_lateral),2)
                n_remove_idx_start(n_group,n_subject,n_lateral) = (n_remove_idx(n_group,n_subject,n_lateral)-1)/2;
                n_remove_idx_end(n_group,n_subject,n_lateral) = (n_remove_idx(n_group,n_subject,n_lateral)+1)/2;
            else
                n_remove_idx_start(n_group,n_subject,n_lateral) = n_remove_idx(n_group,n_subject,n_lateral)/2;
                n_remove_idx_end(n_group,n_subject,n_lateral) = n_remove_idx(n_group,n_subject,n_lateral)/2;
            end
            tmp = [(n_remove_idx_start(n_group,n_subject,n_lateral)+1):(NUM_Depth(n_group,n_subject,n_lateral)-n_remove_idx_end(n_group,n_subject,n_lateral))];
            n_selected_idx(n_group,n_subject,n_lateral,1:N_Selected) = tmp;
            Selected_Depth{n_group,n_subject,n_lateral} = Range_Depth{n_group,n_subject,n_lateral}(tmp(1):tmp(end));
            MIN_Selected_Depth(n_group,n_subject,n_lateral) = min(Selected_Depth{n_group,n_subject,n_lateral});        
            MAX_Selected_Depth(n_group,n_subject,n_lateral) = max(Selected_Depth{n_group,n_subject,n_lateral});
        end % lateral
    end % group
end % subject

%% select subjects
subject_selected_side = ones(2,60,2);
for n_subject=1:max(NewParamArray.subject)
    for n_group=1:2
        for n_lateral=1:2
            n_numeachdepth = Num_EachDepth{n_group,n_subject,n_lateral};
            tmp_idx = squeeze(n_selected_idx(n_group,n_subject,n_lateral,:));
            tmp = n_numeachdepth(tmp_idx);
            if sum(tmp)>N_Selected
                subject_selected_side(n_group,n_subject,n_lateral)=0;
            end
        end
    end
end
subject_selected_L = subject_selected_side(:,:,1);
subject_selected_R = subject_selected_side(:,:,2);
subject_selected = subject_selected_L.*subject_selected_R;
% LA36, GA29

%% select files
for n_file = 1:numel(newfilename)
    file_selected(n_file,1) = 0;
    s_group = NewParamArray.group(n_file);
    if strcmp(s_group,'LA')
        n_group=1;
    elseif strcmp(s_group,'GA')
        n_group=2;
    end
    n_subject = NewParamArray.subject(n_file);
    s_lateral = NewParamArray.lateral(n_file);
    if strcmp(s_lateral,'L')
        n_lateral=1;
    elseif strcmp(s_lateral,'R')
        n_lateral=2;
    end
    s_depth = NewParamArray.depth(n_file);
    % determine whether this file should be selected
    if ismember(s_depth,Selected_Depth{n_group,n_subject,n_lateral})
        file_selected(n_file,1) = 1;
    end
    if subject_selected(n_group,n_subject)==0
        file_selected(n_file,1) = 0;
    end
end

%% copy files
for n_file = 1:numel(newfilename)
    if file_selected(n_file,1)
        source_file = ['MATFilesNew\',newfilename{n_file},'.mat'];
        destination_file = ['MATFilesSelected\',newfilename{n_file},'.mat'];        
        copyfile(source_file,destination_file);
    end
end