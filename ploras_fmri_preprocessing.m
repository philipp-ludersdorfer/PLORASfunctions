function ploras_fmri_preprocessing(subjects,searchdir,TA)
% PLORAS_FMRI_PREPROCESSING(SUBJECTS,SEARCHDIR,TA) is the automatic preprocessing 
% pipeline for the PLORAS 3 fMRI paradigm. Preprocessing steps include realignment 
% & unwarping of the functional MRI images, co-registration between functional 
% and structural MRI images, normalisation to standard MNI space, and smoothing. 
%
% Required inputs: 
%
% SUBJECTS: The IDs of the subjects to include. Either a single subject ID
% (e.g. 'PS0111') or a textfile containing a list of subject IDs (in a single
% column). 
%
% SEARCHDIR: The directory which includes the subject directories. (please 
% provide the full path in inverted commas!).
% 
% TA (optional): Set to 1, this will additionally run a time series analysis 
% before and after the other preprocessing steps. This analysis visualises 
% the mean and the variance of signal intensities per scan. The output of the 
% time series analysis is printed to a .ps file in the subject-specific 
% preprocessed_data folder. 
% 
% (TA requires that the Import_Archive toolbox folder (see subfunctions folder 
% in ploras_functions) is stored in the SPM toolbox directory!)
%
% Advanced options:
%
% To switch off smoothing, make changes to the image resolution or the size 
% of the smoothing kernel it is necessary to make manual edits to this function. 
% To do this type ‘edit ploras_fmri_preprocessing.m’ in the Matlab command window
% and make the desired changes in the now opened script (go to the 'MANUAL
% EDITS' section).
%
% Before running the edited script make sure you save it under a different 
% name (otherwise you will overwrite the default settings)!
% 
% Philipp Ludersdorfer (last edited on 29/11/2016)

%%
%%%%%%%%%%%%%%%%%%%%   M A N U A L   E D I T S   %%%%%%%%%%%%%%%%%%%%%%%%%%

% Smoothing
Smooth = 1; % Set to 0 to omit smoothing step
fwhm = [6 6 6]; % Size of smoothing kernel in mm (default = [6 6 6])

% Image resolution
f_res = [3 3 3]; % Resolution of functional images in mm (default = [3 3 3])
s_res = [1 1 1]; % Resolution of structural images in mm (default = [1 1 1])

% If you are making changes in this section of the script make sure you
% SAVE the script under a different name before running it!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
if(nargin<3) TA=0; end

spm('defaults', 'FMRI') % start SPM
fail = {}; fsi = 1;

%%  Loop over subjects
for crun = 1:length(FinalFolders)
    disp('# # # # # # # # # # # # # # # # # # # # # # # # # #')
    disp(['Running preprocessing for subject ' FinalSubjects{subj}])
    disp('# # # # # # # # # # # # # # # # # # # # # # # # # #')
    
    % Setup input data (= image files)
    subjpath = fullfile(searchdir,FinalFolders{crun},'functional','preprocessed_data');
    [~,dat] = spm_select('List',subjpath); % select fMRI data folders (usually 14 folders for each of the tasks)
    dat = cellstr(dat); % convert to cellstring
    
    % Reorder task folders appropriately
    [~,datsort] = strtok(dat,'.');
    datsort(cellfun(@isempty, datsort)) = [];
    [~,i] = sort(cellfun(@str2num,cellfun(@(x) x(2:end),datsort,'UniformOutput',false))); %you need to change the '2' to '3' if there's a 'S' proceeding the number of the scan run (e.g. S3)
    dat = dat(i);

    % Setup inputs variable with functional and structural image names
    for j=1:length(dat)
        inputs{j,1} = cellstr(spm_select('FPList',fullfile(subjpath,dat{j}),'^f.*\.nii$'));
    end
    inputs{length(dat)+1,1} = cellstr(spm_select('FPList',fullfile(searchdir,FinalFolders{crun},'structural'),'^sM.*\.nii$'));
    
    
    %% SETUP PROCESSING BATCH %%
    step = 1;
    
    % Time Series Analysis (Pre Preprocessing)
    if TA == 1
        matlabbatch{step}.spm.tools.tsdiffana_tools{1}.tsdiffana_timediff.imgs = inputs(1:length(dat))';
        matlabbatch{step}.spm.tools.tsdiffana_tools{1}.tsdiffana_timediff.vf = false;
        for j = 1:length(dat)
            matlabbatch{step+1}.spm.tools.tsdiffana_tools{1}.tsdiffana_tsdiffplot.tdfn(j) = cfg_dep(['Analyse Time Series: Timeseries Analysis Data File (' num2str(j) ')'], substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tdfn', '()',{j}));
        end
        step = step + 2; % update step
    end
 
    % Realign & unwarp
    for j = 1:length(dat)
        matlabbatch{step}.spm.spatial.realignunwarp.data(j).scans = inputs{j};
        matlabbatch{step}.spm.spatial.realignunwarp.data(j).pmscan = '';
    end
    matlabbatch{step}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
    matlabbatch{step}.spm.spatial.realignunwarp.eoptions.sep = 4;
    matlabbatch{step}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
    matlabbatch{step}.spm.spatial.realignunwarp.eoptions.rtm = 0;
    matlabbatch{step}.spm.spatial.realignunwarp.eoptions.einterp = 2;
    matlabbatch{step}.spm.spatial.realignunwarp.eoptions.ewrap = [0 1 0];
    matlabbatch{step}.spm.spatial.realignunwarp.eoptions.weight = '';
    matlabbatch{step}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
    matlabbatch{step}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
    matlabbatch{step}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
    matlabbatch{step}.spm.spatial.realignunwarp.uweoptions.jm = 0;
    matlabbatch{step}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
    matlabbatch{step}.spm.spatial.realignunwarp.uweoptions.sot = [];
    matlabbatch{step}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
    matlabbatch{step}.spm.spatial.realignunwarp.uweoptions.rem = 1;
    matlabbatch{step}.spm.spatial.realignunwarp.uweoptions.noi = 5;
    matlabbatch{step}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
    matlabbatch{step}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
    matlabbatch{step}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
    matlabbatch{step}.spm.spatial.realignunwarp.uwroptions.wrap = [0 1 0];
    matlabbatch{step}.spm.spatial.realignunwarp.uwroptions.mask = 1;
    matlabbatch{step}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
    
    % Co-register
    matlabbatch{step+1}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{step}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
    matlabbatch{step+1}.spm.spatial.coreg.estimate.source = inputs{length(dat)+1}; % Structural Scan Input
    matlabbatch{step+1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{step+1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{step+1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{step+1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{step+1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    % Normalisation 
    matlabbatch{step+2}.spm.spatial.normalise.est.subj.vol(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{step+1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{step+2}.spm.spatial.normalise.est.eoptions.biasreg = 0.0001;
    matlabbatch{step+2}.spm.spatial.normalise.est.eoptions.biasfwhm = 60;
    matlabbatch{step+2}.spm.spatial.normalise.est.eoptions.tpm = {fullfile(spm('Dir'),'tpm','TPM.nii')};
    matlabbatch{step+2}.spm.spatial.normalise.est.eoptions.affreg = 'mni';
    matlabbatch{step+2}.spm.spatial.normalise.est.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{step+2}.spm.spatial.normalise.est.eoptions.fwhm = 0;
    matlabbatch{step+2}.spm.spatial.normalise.est.eoptions.samp = 3;
    matlabbatch{step+3}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Normalise: Estimate: Deformation (Subj 1)', substruct('.','val', '{}',{step+2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','def'));
    matlabbatch{step+3}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{step+1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{step+3}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{step+3}.spm.spatial.normalise.write.woptions.vox = s_res;
    matlabbatch{step+3}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{step+3}.spm.spatial.normalise.write.woptions.prefix = 'w';
    matlabbatch{step+4}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Normalise: Estimate: Deformation (Subj 1)', substruct('.','val', '{}',{step+2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','def'));
    for j = 1:length(dat)
    matlabbatch{step+4}.spm.spatial.normalise.write.subj.resample(j) = cfg_dep(['Realign & Unwarp: Unwarped Images (Sess ' num2str(j) ')'], substruct('.','val', '{}',{step}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{j}, '.','uwrfiles'));
    end
    matlabbatch{step+4}.spm.spatial.normalise.write.subj.resample(length(dat)+1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{step}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
    matlabbatch{step+4}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{step+4}.spm.spatial.normalise.write.woptions.vox = f_res;
    matlabbatch{step+4}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{step+4}.spm.spatial.normalise.write.woptions.prefix = 'w';
    
    % Smoothing
    if Smooth == 1
        matlabbatch{step+5}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{step+4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
        matlabbatch{step+5}.spm.spatial.smooth.fwhm = fwhm;
        matlabbatch{step+5}.spm.spatial.smooth.dtype = 0;
        matlabbatch{step+5}.spm.spatial.smooth.im = 0;
        matlabbatch{step+5}.spm.spatial.smooth.prefix = 's';
    end
    
    % Time Seris Analysis (Post-preprocessing)
    if TA == 1
        if Smooth == 1
            matlabbatch{step+6}.spm.tools.tsdiffana_tools{1}.tsdiffana_timediff.imgs{1}(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{step+5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
            matlabbatch{step+6}.spm.tools.tsdiffana_tools{1}.tsdiffana_timediff.vf = false;
            matlabbatch{step+7}.spm.tools.tsdiffana_tools{1}.tsdiffana_tsdiffplot.tdfn(1) = cfg_dep('Analyse Time Series: Timeseries Analysis Data File (1)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tdfn', '()',{1}));
            matlabbatch{step+7}.spm.tools.tsdiffana_tools{1}.tsdiffana_tsdiffplot.doprint = true;
        else
            matlabbatch{step+5}.spm.tools.tsdiffana_tools{1}.tsdiffana_timediff.imgs{1}(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{step+4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
            matlabbatch{step+5}.spm.tools.tsdiffana_tools{1}.tsdiffana_timediff.vf = false;
            matlabbatch{step+6}.spm.tools.tsdiffana_tools{1}.tsdiffana_tsdiffplot.tdfn(1) = cfg_dep('Analyse Time Series: Timeseries Analysis Data File (1)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tdfn', '()',{1}));
            matlabbatch{step+6}.spm.tools.tsdiffana_tools{1}.tsdiffana_tsdiffplot.doprint = true;
        end
    end
    
    %% RUN PREPROCESSING BATCH %%
    cd(subjpath);
    save(['P3_prepro_job_' datestr(date,'yyyymmmdd')],'matlabbatch') % save job as .mat file (can be later opened in SPM batch editor)
    
    try % Try to run batch
        spm_jobman('run', matlabbatch); % run jobman
    catch % If error occurs
        fail{fsi} = [FinalSubjects{crun} ': SPM batch failed'];
        fsi = fsi + 1;
    end
    
    try % Try to rename output .ps file (only exists if TS analysis was performed!)
        log = spm_select('FPList',subjpath,['^spm_' datestr(date,'yyyymmmdd') '.ps']);
        movefile(log,[subjpath '\P3_prepro_plots_' datestr(date,'yyyymmmdd') '.ps']) 
    catch
    end
    
    clear matlabbatch inputs step log
end

%% Summary
disp('# # # # # # # # # # # # # # # # # # # # # # # # # #')
disp('Congratulations, the preprocessing is done!')
disp('# # # # # # # # # # # # # # # # # # # # # # # # # #')
if fsi ~= 1 % Errors
    disp('')
    disp('The following errors occurred:')
    for i = 1:length(fail)
        disp(fail{i});
    end
end


