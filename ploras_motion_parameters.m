function [motionmat motvec] = ploras_motion_parameters(subjects,tasks,searchdir,output)
% PLORAS_MOTION_PARAMETERS(SUBJECTS,TASKS,SEARCHDIR,OUTPUT) extracts task-
% specific head motion parameters (path lengths) for the PLORAS 3 fMRI paradigm
% and saves them in an Excel spreadsheet.
%
% Required inputs:
% 
% SUBJECTS: either a single subject ID in inverted commas (e.g. 'PS0111') or
% a text file containing a list of subject IDs (in a single column). If the same
% subject had more than one scan, the script will automatically detect
% this. This means that the ID does not have to be included more than once!
% 
% TASKS: a numerical vector indicating the indices of the tasks of the P3 fMRI 
% paradigm to be included. Examples are: 1:14 (for all tasks), [1 8 5 6] (for 
% a subset of the tasks in a particular order).
%
% (HINT: If you don’t know the index/indices of your task/s of interest, type 
% ploras_fmri_conditions in the MATLAB command line. This will list all tasks 
% and their corresponding indices.)
%
% SEARCHDIR: The directory in which the subject directories can be found
% (please provide the full path in inverted commas!).
%
% OUTPUT: Name of the output Excel spreadsheet, e.g. '30_patients_summary'.
% If you don't provide a full path, the file will be stored in the current
% working directory.
%
% Philipp Ludersdorfer (last modified 02/12/2016)
% Based on previous versions by Mohamed Seghier (2007), Suz Prejawa and
% 'Oiwi Parker Jones (2014)


%% Check input arguments
% subject list and searchdir
if ischar(subjects) && ischar(searchdir) % if string
    if strcmp(subjects(end-2:end),'txt')
        try
            SubjList = textread(subjects,'%s');
        catch
            error('Subject list text file is not an existing filename')
        end
    else
        SubjList{1} = subjects;
    end
    subfolds = dir(searchdir);
    AllFolders = {subfolds([subfolds.isdir]).name};
    fin = 1;
    for i = 1 : length(SubjList)
        dirind = find(~cellfun(@isempty,strfind(AllFolders,SubjList{i})));
        for j = 1 : length(dirind)
            FinalFolders(fin) = AllFolders(dirind(j));
            FinalSubjects(fin) = SubjList(i);
            fin = fin + 1;
        end
    end
end
% tasks
if(isnumeric(tasks))
    TaskList = ploras_fmri_conditions(tasks);
else
    error('Invalid task specfication!')
end
% output
if(~ischar(output))
    error('Invalid output filename!')
end

motionmat = zeros(length(FinalSubjects),length(tasks)*6); % initialise individual motion parameter vector
v0 = 2*[10;10;10] ; % reference point for rotation ISD calculation (see Yoo et al. 2005, Neuroscience Res)
fail = {}; fsi = 1;

% Loop over subjects
disp('Running motion parameter extraction ...')
for subj=1:numel(FinalFolders)
    
    disp(['Subject: ' FinalSubjects{subj}])
   
    % fMRI data folders
    [~,dat] = spm_select('List',fullfile(searchdir,FinalFolders{subj},'functional','preprocessed_data'));
    dat = cellstr(dat);
    % and reorder them appropriately
    [~,datsort] = strtok(dat,'.');
    [~,i] = sort(cellfun(@str2num,cellfun(@(x) x(2:end),datsort,'uniformoutput',false)));
    dat = dat(i);
    
    if numel(dat) < numel(TaskList) % Check if there are enought fMRI data folders
        fail{fsi} = [FinalSubjects{subj} ': There are less functional folders than specified tasks'];
        fsi = fsi + 1;
        continue 
    end
    
    % loop over tasks/fmri sessions
    for s = 1 : numel(tasks)
        
        % Check and open session-specific motion parameter file
        rpfile = spm_select('FPList',fullfile(searchdir,FinalFolders{subj},'functional','preprocessed_data',dat{tasks(s)}),'^rp.*.txt');
        if ~exist(rpfile,'file')
            fail{fsi} = [FinalSubjects{subj} ': No existing head motion file'];
            fsi = fsi + 1;
            continue
        end
        rp = load(rpfile);
        
        % Remove drift from input signal
        rp = moh_rmdrift(rp) ;
        
        %% Extract motion parameters
        % 1. Translation (first 3 columns of rp file)
        % Calculate interscan displacement (ISD) in mm (= Euclidean distance 
        % between two successive scans (modified from Zou etal. 2005 Neuroimage)
        T_ISD = sqrt( (diff(rp(:,1)).^2) + (diff(rp(:,2)).^2) +...
            (diff(rp(:,3)).^2) ) ; 
        T_mean(s) = mean(T_ISD); % Mean ISD (path length) in mm
        T_max(s) = max(T_ISD); % Maximal ISD 
        T_sum(s) = sum(T_ISD); % Total ISD
        
        % 2. Rotation (column 4 to 6 of rp file)
        % Convert degrees to mm and calculate interscan displacement (ISD) 
        % (for rationale see Yoo et al. 2005 Neurosci Res)
        for j=1:length(rp(:,1))
            ax = rp(j,4) ;
            ay = rp(j,5) ;
            az = rp(j,6) ;
            Rx = [1, 0, 0; 0, cos(ax), sin(ax); 0, -sin(ax), cos(ax)] ;
            Ry = [cos(ay), 0, sin(ay); 0, 1, 0; -sin(ay), 0, cos(ay)] ;
            Rz = [cos(az), sin(az), 0; -sin(az), cos(az), 0; 0 , 0,1] ;
            Rxyz(:,:,j) = Rx * Ry * Rz ;
            Rotat(j,:)  = (Rxyz(:,:,j) * v0)';
        end
        R_ISD = sqrt( (diff(Rotat(:,1)).^2) + (diff(Rotat(:,2)).^2) +...
            (diff(Rotat(:,3)).^2) ) ;
        R_mean(s) = mean(R_ISD); % Mean ISD (path length) in mm
        R_max(s) = max(R_ISD); % Maximal ISD 
        R_sum(s) = sum(R_ISD); % Total ISD
        motvec{subj,s} = {T_ISD R_ISD};
    end
    % Save parameters into results matrix
    motionmat(subj,:) = [T_mean, R_mean, T_sum, R_sum, T_max, R_max]; 
end

%% Write Excel spreadsheet
% Setup cell array
cols = {'Trans Mean ISD','Rot Mean ISD','Trans Total ISD','Rot Total ISD','Trans Max ISD','Rot Max ISD'};
repcols = repmat(cols,length(TaskList),1); 
cols = repcols(:)';
MotionArray(1,:) = {'', repcols{1:end}} ;
taskcols = repmat(TaskList',1,6);
MotionArray(2,:) = {'Subject',taskcols{1:end}};
MotionArray(3:(2+length(FinalSubjects)),1) = FinalSubjects;
MotionArray(3:(2+length(FinalSubjects)),2:(size(motionmat,2)+1)) = num2cell(motionmat);
% Write to file
xlswrite([output '.xlsx'], MotionArray);

%% Display warnings
disp('')
disp('Motion parameter extraction is done!')
if fsi ~= 1 % Errors
    disp('')
    disp('The following errors occurred:')
    for i = 1:length(fail)
        disp(fail{i});
    end
end

end


function evolF = moh_rmdrift(evol)
% remove a drift from an input signal

if size(evol,1) < size(evol,2), evol = evol'; end

NVolume = size(evol, 1);
K = 5 ; % Dimension de l'espace vectoriel o.n dans lequel sera projeté le signal d'évolution temporelle d'un voxel
CosBase = zeros (NVolume, K); % composantes en COS (prises en compte par SPM)
Theta = zeros (K, 1); % coefficient de chacun des COS
Alpha = ((0:NVolume-1) / (NVolume))';
for r = 1:K
    Phi (:, r) = cos((r) * pi * Alpha);
end;
PInvPhi = (Phi'*Phi)^(-1);

for i=1:size(evol,2)
    % Correction de la derive de ligne de base des evolutions des pixels
    %-------------------------------------------------------------------
    Theta = PInvPhi * Phi' * (evol(:,i) - mean(evol(:,i)));
    evolF(:,i) = evol(:,i) - Phi*Theta; % Evolution temporelle filtrée (- la ligne de base)
end
end