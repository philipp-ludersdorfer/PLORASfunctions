function ploras_data_import(src,dest,tweak)
% PLORAS_DATA_IMPORT(SRC,DEST,TWEAK) imports and extracts archives from 
% specified source directories on 'Charm' into specified destination directories 
% and converts images into '.nii' format. The function additionally sets up a 
% folder structure in destination directory for functional and structural 
% images and moves image files into their respective folders. 
%
% Can be used in two ways:
%
% 1. When processing only one subject: Two strings are required as input. One string 
% is the source directory on 'Charm' and the other the destination directory (both 
% have to include the full path).
% e.g. ploras_data_import('W:\trio\2013\Example','Z:\PLORAS3_fMRI_01\Example')
%
% 2. When processing more than one subject: A single textfile is required as input. 
% The textfile has to include two columns. The first column should contain the source
% directory names on 'Charm' and the second column should contain the corresponding 
% destination directories (both have to include the full path).
% e.g. ploras_data_import('folders.txt')
%
% Requirement: 
% Each destination directory must contain a run_info.txt file containing the run 
% number - runtype correspondences separated by line breaks
% (e.g. 
% 02 Fieldmap
% 03 Functional
% 04 Resting_state
% 05 Structural
% ...)
%
% TWEAK (optional): Setting the third input argument to 1 allows to import data 
% from a different scan as the "main" data of a subject (i.e. source and 
% destination do not have to include the same scanID).
%
% Philipp Ludersdorfer (last modified 28/11/2016)

if ~exist('dest') || isnumeric(dest) % 1 input argument -> input has to be an existing text file
        if ~exist(src,'file') 
            disp(['ERROR: ' src ' is not an existing text file!']);
            return
        else
            [source,destination] = textread(src,'%s %s'); % read in text file
            if isnumeric(dest)
                tweak = dest;
            end
        end
else
        source{1} = src; % source directory 
        destination{1} = dest; % destination directory
end

for i = 1:length(source) % Check if source and destination only include existing directories
   if ~exist(source{i},'dir') % In case source{i} is not an existing directory
       disp(['ERROR: ' source{i} ' is not an existing folder!'])
       return
   end
   if ~exist(destination{i},'dir') % In case destination{i} is not an existing directory
       disp(['ERROR: ' destination{i} ' is not an existing folder!'])
       return
   end      
end

% If foldernames don't end with '\' -> append '\'
for i = 1:length(source)
    if ~strcmp(source{i}(end),'\')
        source{i}(end+1)='\';
    end
    if ~strcmp(destination{i}(end),'\')
        destination{i}(end+1)='\';
    end
    id.dest{i} = destination{i}(end-4:end-1); % save scan id from destination path (used for comparison with source scan id below)
end
clear src dest

for i = 1:length(source) % loop across inputs
    disp(['## Processing subject ' num2str(i) '/' num2str(length(source)) ' ##'])
    tic; % start clock 
    mkdir(destination{i},'functional\raw_data'); % create Functional\raw_data folder in destination folder 
    
    %% IMPORT AND CONVERSION
    if ~exist([destination{i} 'run_info.txt'],'file') % Check if run_info.txt exits
        disp(['ERROR: No run_info.txt file found in ' destination{i}]);
        return
    else
        [ses.number, ses.name] = textread([destination{i} 'run_info.txt'],'%u %s'); % Open run_info.txt in DESTFOLDER
    end
    sourcefiles = dir([source{i} '*.tar']); % sourcefiles = all .tar files in source folder
       
    for j = 1 : length(sourcefiles) % loop over .tar files
        dotid = strfind(sourcefiles(j).name,'.'); % look for '.' in filenames
        if j == 1 % For first sourcefile check if scan id matches with scan id of destination foldername
            usid = strfind(sourcefiles(j).name,'_'); % look for '_' in filenames
            if ~isempty(usid) 
                id.source{i} = sourcefiles(j).name(1:usid(1)-1); % extract scan id from source filename
            else % in case there is no '_' in source filename
                id.source{i} = sourcefiles(j).name(1:dotid(1)-1); % extract scan id from source filename
            end
            if ~isequal(id.source{i}(end-3:end),id.dest{i}) % In case scan ids do not match between source and destination
                if(tweak)
                    disp(['Scan ID does not match between source and destination folder for ' id.source{i}])
                    disp('I will continue though, because I simply do not care!')
                else
                    disp(['ERROR: Scan ID does not match between source and destination folder for ' id.source{i}])
                    return
                end
            end
        end % end: check if scan ids match
        if strcmp(sourcefiles(j).name(dotid(1)+1),'S') % check if session number starts with 'S'
            scannum = sourcefiles(j).name(dotid(1)+2:dotid(2)-1); % set session number
        else
            scannum = sourcefiles(j).name(dotid(1)+1:dotid(2)-1); % see session number (when it doesn't start with 'S')
        end
        if ismember(str2double(scannum),ses.number) % Check if current session number is of interest (= compare with numbers from run_info.txt)
            copyfile([source{i} sourcefiles(j).name],[destination{i} 'functional\raw_data']) % first copy .tar file from charm (otherwise script is likely to crash due to network connection failure)
            run .\subfunctions\Import_Archive([destination{i} 'functional\raw_data\' sourcefiles(j).name],[destination{i} 'functional\raw_data']); % Import and convert files into 'raw_data' folder in DESTFOLDER
            delete([destination{i} 'functional\raw_data\' sourcefiles(j).name]); % delete copied .tar file
            ses.folder{ses.number==str2double(scannum)} = [destination{i} 'functional\raw_data\' sourcefiles(j).name(1:end-4) '\']; % save outputfolder for later
        end
    end
    fprintf('Import and conversion are done!\n')
    
    %% Move Functional data, Fieldmap(s), Structural data, and Resting state data
    mkdir([destination{i} 'functional\dummies']);
    for j = 1 : length(ses.name) % loop over sessions
        if strcmp(ses.name{j},'Functional') % if session is functional
            if ses.number(j) < 10 % if session id is 1-9 prepend folder name with 0 
                foldername = ([id.source{i} '.0' num2str(ses.number(j))]);
            else
                foldername = ([id.source{i} '.' num2str(ses.number(j))]);
            end
            mkdir([destination{i} 'functional\preprocessed_data'], foldername); % create folder with new name
            filelist = spm_select('FPList',ses.folder{j},'.nii'); % Get a list of all .nii files in session-specific folder in raw_data
            imgfiles = cellstr(filelist); % get filelist in a useful format (cell instead of character matrix)
            for imgcount = 1 : length(imgfiles) % loop over image (=.nii) files
                if imgcount < 6 % images 1-5
                    movefile(imgfiles{imgcount},[destination{i} 'functional\dummies']); % copy into 'Functional\Dummies'
                else % images 6-end
                    copyfile(imgfiles{imgcount},[destination{i} 'functional\preprocessed_data\' foldername]); % copy into 'Functional\Functional_Data_Preprocessing\FunctionalX'
                end
            end
            movefile(ses.folder{j},[destination{i} 'functional\raw_data\' foldername]); % change name of raw_data folder (Append 0)!
        elseif strcmp(ses.name{j},'Fieldmap') % if session is a fieldmap
            movefile([ses.folder{j} 'sM*'],[destination{i} 'functional\fieldmap']); % move files to Fieldmap directory
            rmdir(ses.folder{j}) % remove raw data directory
        elseif strcmp(ses.name{j},'Fieldmap2') % if there was a second fieldmap
            movefile([ses.folder{j} 'sM*'],[destination{i} 'functional\fieldmap2']);
            rmdir(ses.folder{j}) % remove raw data directory
        elseif strcmp(ses.name{j},'Fieldmap3') % if there was a third fieldmap
            movefile([ses.folder{j} 'sM*'],[destination{i} 'functional\fieldmap3']);
            rmdir(ses.folder{j}) % remove raw data directory
        elseif strcmp(ses.name{j},'Structural') % if session is structural
            movefile(ses.folder{j},[destination{i} 'structural']);
        elseif strcmp(ses.name{j},'Resting_state') % if session is resting state
            movefile(ses.folder{j},[destination{i} 'functional\raw_data\' 'resting_state']); % Rename folder into "Resting_state"
            filelist = spm_select('FPList',[destination{i} 'functional\raw_data\' 'resting_state'],'.nii'); % Get a list of all .nii files in session-specific folder in raw_data
            imgfiles = cellstr(filelist); % get filelist in a useful format (cell instead of character matrix)
            for imgcount = 1:5 % 5 Dummy scans
                movefile(imgfiles{imgcount},[destination{i} 'functional\dummies']); % copy Dummies into 'Functional\Dummies'
            end
        end
    end
    clear dotid filelist foldername imgcount imgfiles ses scannum sourcefiles usid % just to be sure
    fprintf('Setup of folder structure and moving of data done!\n')
    toc
    fprintf('\n')
end
fprintf('## Well done! Everything is in its right place! ##\n')
