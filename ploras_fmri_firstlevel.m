function ploras_fmri_firstlevel(subjects,tasks,searchdir,complex)
% PLORAS_FMRI_FIRSTLEVEL(SUBJECTS,TASKS,SEARCHDIR,COMPLEX) runs the first 
% level (i.e. subject-specific) analysis for the PLORAS 3 fMRI paradigm. 
% The resulting SPM.mat is stored in a newly created 'first_level' folder 
% in the subject directory.
%
% Required inputs:
% 
% SUBJECTS: either a single subject ID in inverted commas (e.g. 'PS0111’) or 
% a text file containing a list of subject IDs (in a single column). 
%
% TASKS: a numerical vector indicating the indices of the tasks of the P3 
% fMRI paradigm to be included. Examples are: 1:14 for all tasks, 1:5 or 
% [1 2 3 4 5] for the first 5 tasks, or [1:3 5:14] when task no. 4 is missing. 
%
% (Hint: If you don’t know the index/indices of your task/s of interest, type 
% ploras_fmri_conditions in the MATLAB command window. This will list all tasks 
% and their corresponding indices.)
%
% SEARCHDIR: The directory in which the subject directories can be found
% (please provide the full path in inverted commas!).
%
% COMPLEX: the analysis is either based on the simple model (0 = default) 
% or the complex model (1).
% 
%   Simple model:
%   - 4 regressors per task: instruction, correct, incorrect, missing
%	- 1 contrast per task: correct > fixation baseline
%   - Results are stored in ‘\first_level\simple_model’
%
%	Complex model:
%	- 7 regressors per task: instruction, 4 x correct, incorrect, missing
%	- 5 contrasts per task: all correct > fixation baseline + each condition > fixation baseline
%	- Results are stored in ‘\first_level\complex_model’
%
% Philipp Ludersdorfer (last modified 01/12/2016)


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
    [TaskList CondList] = ploras_fmri_conditions(tasks,1);
else
    error('Invalid task specfication!')
end
% complex
if nargin < 4; complex = 0; end
if complex == 0
    statdir = 'simple_model';
else
    statdir = 'complex_model';
end

spm('defaults','fmri');
spm_jobman('initcfg');

fail = {}; fsi = 1; wrngs = {}; wrn = 1;

for subj=1:numel(FinalFolders)
    disp('# # # # # # # # # # # # # # # # # # # # # # # # # #')
    disp(['Running first-level analysis for subject ' FinalSubjects{subj}])
    disp('# # # # # # # # # # # # # # # # # # # # # # # # # #')
    
    % Open behavioural XL file
    xlsfile = fullfile(searchdir,FinalFolders{subj},'behavioural',[FinalFolders{subj} '.xlsx']);
    if ~exist(xlsfile,'file')
        fail{fsi} = [FinalSubjects{subj} ': No existing behavioural XL file'];
        fsi = fsi + 1;
        continue 
    end
    
    %% 1. MODEL SPECIFICATION
    clear matlabbatch
    mkdir(fullfile(searchdir,FinalFolders{subj},'first_level'),statdir);
    matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(fullfile(searchdir,FinalFolders{subj},'first_level',statdir));
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 3.08;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 44;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 22;
    
    % fMRI data folders
    [~,dat] = spm_select('List',fullfile(searchdir,FinalFolders{subj},'functional','preprocessed_data'));
    dat = cellstr(dat);
    % and reorder them appropriately
    [~,datsort] = strtok(dat,'.');
    [~,i] = sort(cellfun(@str2num,cellfun(@(x) x(2:end),datsort,'uniformoutput',false)));
    dat = dat(i);
    
    if numel(dat) < numel(TaskList)
        fail{fsi} = [FinalSubjects{subj} ': There are less functional folders than specified tasks'];
        fsi = fsi + 1;
        continue 
    end
    
    % get XLS sheets name
    [~,sheets] = xlsfinfo(xlsfile);
    
    for s = 1 : numel(tasks)
        
        % FMRI files
        matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans = cellstr(spm_select('FPList',fullfile(searchdir,FinalFolders{subj},'functional','preprocessed_data',dat{s}),'^swufM.*\.nii$'));
        
        % XLS data
        [~, ~, raw] = xlsread(xlsfile,sheets{tasks(s)});
        raw = raw(2:end,[3:5 8:9]);
        raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
        R = cellfun(@(x) ~isnumeric(x),raw);
        raw(R) = {NaN};
        oart = cell2mat(raw);
        oart(find(all(isnan(oart),2)),:) = [];
        
        % convert new scoring system into old one for task 2 and 4
        if tasks(s) == 2 || tasks(s) == 4
            oart(oart(:,3) == 1 | oart(:,3) == 2 | oart(:,3) == 3,3) = 1;
            oart(oart(:,3) == 4 | oart(:,3) == 5,3) = 2;
        end
        
        % EVENTS ...
        % 1. Instruction
        matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(1).name = 'instruction';
        matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(1).onset = oart(isnan(oart(:,3)),1) / 1000;
        if isempty(matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).onset)
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(1).onset = numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08;
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(1).duration = 0;
        % 2. Correct trials
        if ~complex
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end+1).name = 'correct';
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = oart(oart(:,3) == 2,1) / 1000;
            if isempty(matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset)
                matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08;
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).duration = 0;
        elseif complex
            oart2 = oart(oart(:,3) == 2,:);
            onsets = oart2(:,1) / 1000;
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end+1).name = CondList{s}{1};
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = onsets(oart2(:,4) == 1 & oart2(:,5) == 1);
            if isempty(matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset)
                matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08;
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).duration = 0;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end+1).name = CondList{s}{2};
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = onsets(oart2(:,4) == 2 & oart2(:,5) == 1);
            if isempty(matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset)
                matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08;
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).duration = 0;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end+1).name = CondList{s}{3};
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = onsets(oart2(:,4) == 1 & oart2(:,5) == 2);
            if isempty(matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset)
                matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08;
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).duration = 0;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end+1).name = CondList{s}{4};
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = onsets(oart2(:,4) == 2 & oart2(:,5) == 2);
            if isempty(matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset)
                matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08;
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).duration = 0;
        end
        % 3. Incorrect (or self-corrected) trials
        matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end+1).name = 'incorrect';
        matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = oart(oart(:,3) == 1 | oart(:,3) == 1.5,1) / 1000; % beware
        if isempty(matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset)
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08;
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).duration = 0;
        % 4. No answer
        matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end+1).name = 'no_answer';
        matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = oart(oart(:,3) == 0,1) / 1000;
        if isempty(matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset)
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).onset = numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08;
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(end).duration = 0;
       % %%%% MOVEMENT PARAMETERS %%%%
       % matlabbatch{1}.spm.stats.fmri_spec.sess(s).multi_reg = cellstr(spm_select('FPList',fullfile(searchdir,FinalFolders{subj},'functional','preprocessed_data',dat{s}),'^MovPar*'));
        matlabbatch{1}.spm.stats.fmri_spec.sess(s).hpf = 128;
    end
    
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.mask = {'Z:\PLORAS_fMRI_P3_01_SPM12_FEB2016\Functional\whole_brain_mask.nii'};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    
    %% 2. MODEL ESTIMATION
    matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(searchdir,FinalFolders{subj},'first_level',statdir,'SPM.mat'));
    
    
    %% 3. CONTRAST MANAGER
    matlabbatch{3}.spm.stats.con.spmmat = cellstr(fullfile(searchdir,FinalFolders{subj},'first_level',statdir,'SPM.mat'));
    matlabbatch{3}.spm.stats.con.delete = 1; % Delete existing contrasts? (1 = yes, 0 = no)
    
    if ~complex
        % 3.1. F-map correct
        matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'F-map correct';
        matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec = {kron(eye(length(TaskList)),[0 1 0 0])};
        
        % detect empty 'correct' columns (ie 3rd column in each task) and update 'F-map-correct' contrast matrix
        for s=1:numel(tasks)
            ncond = 2;
            if numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(ncond).onset) == 1 && matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(ncond).onset == numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08
                matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec{1}(:,(s-1)*4+2) = 0;
            end
        end
        matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec{1}(~any(matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec{1},2),:) = [];
        
        % 3.2. F-map incorrect
        matlabbatch{3}.spm.stats.con.consess{end+1}.fcon.name = 'F-map incorrect';
        matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec = {kron(eye(length(TaskList)),[0 0 1 0])};
        
        % detect empty 'incorrect' columns (ie 3rd column in each task) and update 'F-map-incorrect' contrast matrix
        for s=1:numel(tasks)
            ncond = 3;
            if numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(ncond).onset) == 1 && matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(ncond).onset == numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08
                matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec{1}(:,(s-1)*4+3) = 0;
            end
        end
        matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec{end}(~any(matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec{1},2),:) = [];
        
        % 3.3. Baseline T contrasts (i.e. Task [correct trials only!] > Rest)
        for conind = 1:numel(TaskList)
            contask = zeros(1,length(TaskList)); contask(conind) = 1;
            % if there are no correct trials
            if matlabbatch{1}.spm.stats.fmri_spec.sess(conind).cond(2).onset == (numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08)
                matlabbatch{3}.spm.stats.con.consess{end+1}.tcon.name = [TaskList{conind} '_alltrials'];
                wrngs{wrn} = [FinalSubjects{subj} ': Contrast ' num2str(conind+2) ' is based on all trials as there were no correct responses for ' TaskList{conind}];
                wrn = wrn + 1;
                % check if there are 'incorrect' trials or if all responses are missing
                if matlabbatch{1}.spm.stats.fmri_spec.sess(conind).cond(3).onset == (numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08)
                    matlabbatch{3}.spm.stats.con.consess{end}.tcon.convec = kron(contask,[0 0 0 1]);
                elseif matlabbatch{1}.spm.stats.fmri_spec.sess(conind).cond(4).onset == (numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08)
                    matlabbatch{3}.spm.stats.con.consess{end}.tcon.convec = kron(contask,[0 0 1 0]);
                else
                    matlabbatch{3}.spm.stats.con.consess{end}.tcon.convec = kron(contask,[0 0 1 1]);
                end
            % if there are correct trials
            else
                matlabbatch{3}.spm.stats.con.consess{end+1}.tcon.name = TaskList{conind};
                matlabbatch{3}.spm.stats.con.consess{end}.tcon.convec = kron(contask,[0 1 0 0]);
            end
        end
        
    elseif complex
        % 3.1. F-map correct tasks
        matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'F-map-correct Tasks';
        matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec = {kron(eye(length(TaskList)),[0 1 1 1 1 0 0])};
        
        % 3.2. F-map correct conditions
        matlabbatch{3}.spm.stats.con.consess{end+1}.fcon.name = 'F-map-correct Conditions';
        % setup contrast vector only including regressors of interest (no instructions, incorrect, or misses)
        a = eye(length(TaskList)*7); cor = [];
        for i = 1 : length(TaskList)
            cor = [cor; a(((i*7-5):(i*7-2)),:)];
        end
        matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec = {cor}';
        
        % detect empty 'correct' columns in the F-map contrast matrices
        for s=1:numel(tasks)
            for ncond = 2:5
                if numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(ncond).onset) == 1 && matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond(ncond).onset == numel(matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans)*3.08
                    matlabbatch{3}.spm.stats.con.consess{end-1}.fcon.convec{1}(:,(s-1)*4+ncond) = 0;
                    matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec{1}(:,(s-1)*4+ncond) = 0;
                end
            end
        end
        matlabbatch{3}.spm.stats.con.consess{end-1}.fcon.convec{1}(~any(matlabbatch{3}.spm.stats.con.consess{end-1}.fcon.convec{1},2),:) = [];
        matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec{1}(~any(matlabbatch{3}.spm.stats.con.consess{end}.fcon.convec{1},2),:) = [];
        
        % 3.3. Task > Rest t contrasts (averaged across within-task conditions)
        for conind = 1:numel(TaskList)
            contask = zeros(1,length(TaskList)); contask(conind) = 0.25;
            matlabbatch{3}.spm.stats.con.consess{end+1}.tcon.name = TaskList{conind};
            matlabbatch{3}.spm.stats.con.consess{end}.tcon.convec = kron(contask,[0 1 1 1 1 0 0]);
        end
        
        % 3.4. Condition > Rest t contrasts
        for conind = 1:numel(TaskList)
            contask = zeros(1,length(TaskList)); contask(conind) = 1;
            % Condittion 1
            matlabbatch{3}.spm.stats.con.consess{end+1}.tcon.name = CondList{conind}{1};
            matlabbatch{3}.spm.stats.con.consess{end}.tcon.convec = kron(contask,[0 1 0 0 0 0 0]);
            % Condition 2
            matlabbatch{3}.spm.stats.con.consess{end+1}.tcon.name = CondList{conind}{2};
            matlabbatch{3}.spm.stats.con.consess{end}.tcon.convec = kron(contask,[0 0 1 0 0 0 0]);
            % Condition 3
            matlabbatch{3}.spm.stats.con.consess{end+1}.tcon.name = CondList{conind}{3};
            matlabbatch{3}.spm.stats.con.consess{end}.tcon.convec = kron(contask,[0 0 0 1 0 0 0]);
            % Condition 4
            matlabbatch{3}.spm.stats.con.consess{end+1}.tcon.name = CondList{conind}{4};
            matlabbatch{3}.spm.stats.con.consess{end}.tcon.convec = kron(contask,[0 0 0 0 1 0 0]);
        end
    end
    
    %% 4. SAVE JOB AS BATCH AND RUN
    % Save job as batch file....
    save(fullfile(searchdir,FinalFolders{subj},'first_level',statdir,['firstlevel_' statdir '_batch.mat']),'matlabbatch');
    
    % Run job
    try % Try to run batch
        spm_jobman('run', matlabbatch); % run jobman
    catch % If error occurs
        fail{fsi} = ['The SPM batch failed for ' FinalSubjects{subj}];
        fsi = fsi + 1;
    end
    
end

%% Summary
disp('# # # # # # # # # # # # # # # # # # # # # # # # # #')
disp('Congratulations, the first-level analysis is done!')
disp('# # # # # # # # # # # # # # # # # # # # # # # # # #')
if fsi ~= 1 % Errors
    disp('')
    disp('The following errors occurred:')
    for i = 1:length(fail)
        disp(fail{i});
    end
end
if wrn ~= 1
    disp('')
    disp('The following warnings occurred:')
    for i = 1:length(wrngs)
        disp(wrngs{i})
    end
end