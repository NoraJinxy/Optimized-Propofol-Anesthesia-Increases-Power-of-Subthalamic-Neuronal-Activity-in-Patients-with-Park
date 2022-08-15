% this file is used to read text files, extract info from files/filenames, 
% and to write data to MAT files
clear all

%% List all text files
dirnames = {'A_STN';'B1_STN';'B2_STN'};
filenames = [];
for n_dir=1:numel(dirnames)
    filelist = dir(['TXTFiles\',dirnames{n_dir}]);
    for n_file=3:numel(filelist)
        tmp = filelist(n_file).name;
        if strcmp(tmp(end-3:end),'.txt')
            filenames{n_dir,1}{n_file-2,1} = tmp(1:end-4);
        end
    end
end

%% Read text and write to mat
for n_dir=1:numel(dirnames)
    n_dir
    for n_file=1:numel(filenames{n_dir})
        n_file
        txtfilename = ['TXTFiles\',dirnames{n_dir},'\',filenames{n_dir}{n_file},'.txt'];
        data = dlmread(txtfilename);
        x.signal = data(:,2); % MER signals
        x.time = data(:,1); % time index
        x.other = data(:,3:end); % ?
        x.fs = 1000; % sampling rate
        % extract info from file name
        strinfo = filenames{n_dir}{n_file};
        count_ = 0; % find number of "_" in the file name
        for i=1:numel(strinfo)
            if strcmp(strinfo(i),'_')
                count_ = count_+1;
                idx_(count_) = i;
            end
        end
        if count_==5
            x.param1 = dirnames{n_dir}(1:end-4); % A/B?
            x.param2 = strinfo(1:idx_(1)-1); % subject/dir
            x.param3 = strinfo(idx_(1)+1:idx_(2)-1); % left or right
            x.param4 = strinfo(idx_(2)+1:idx_(3)-1); % ?
            x.param5 = strinfo(idx_(3)+1:idx_(4)-1); % ?
            x.param6 = strinfo(idx_(4)+1:idx_(5)-1); % ?
            x.param7 = strinfo(idx_(5)+1:end); % depth
        else
            disp('WRONG FILE NAME')
            STOP
        end

        %% write
        matfilename = ['MATFiles\',dirnames{n_dir}(1:end-4),'_',filenames{n_dir}{n_file},'.mat'];
        % if the character is '-' or "+', it must be replaced by 'n' or 'p' for MATLAB
        % 'n' for negative, 'p' for positive
%         for i=1:numel(matfilename)
%             if strcmp(matfilename(i),'-')
%                 matfilename(i) = 'n';
%             elseif strcmp(matfilename(i),'+')
%                 matfilename(i) = 'p';
%             end
%         end
        save (matfilename,'x');
    end
end
