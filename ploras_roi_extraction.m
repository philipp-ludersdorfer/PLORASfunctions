function roidata = ploras_roi_extraction(subjects,tasks,searchdir,region,stat)
% PLORAS_ROI_EXTRACTION(SUBJECTS,TASKS,REGION,STAT) extracts the task versus
% fixation baseline parameter estimates from voxels within a specified region 
% of interest and returns a summary statistic for each specified subject.
% The resulting variable ROIDATA is a matrix where each subject corresponds 
% to a row and each task corresponds to a column. 
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
% SEARCHDIR: The directory which contains the subject directories can be found 
% (please provide the full path in inverted commas!).
%
% REGION: region of interest. Can be either the name of an existing region of 
% interest file (in inverted commas) or a numerical vector of length 4 specifying
% a spherical region (i.e. [x, y, z, radius]).
%
% STAT: summary statistic calculated across voxels included in the region
% of interest. Default is the first Eigenvariate of a single value
% decomposition algorithm (as standard in SPM). Can be changed to @mean,
% etc. if necessary.
%
% Philipp Ludersdorfer (last updated 05/12/2016)

% subject list and searchdir
if ~ischar(subjects)
    error('Invalid subject specification!')
end
% tasks
if ~isnumeric(tasks)
    error('Invalid task specfication!')
end
% searchdir
if ~ischar(searchdir)
    error('Invalid searchdir specification')
end
% region
if ischar(region)
    if ~exist(region,'file')
        error('Region is not an exiting filename!')
    end
elseif isnumeric(region)
    if length(region)~=4
        error('Invalid region definition! Must be a vector of length 4 [x, y, z, radius])')
    end
else
    error('Invalid region definition! Must be either an existing filename or a coordinates vector!')
end
% stat
if nargin < 5
    stat = @firsteigenvariate;
end

%% Data extraction
for i = 1:length(tasks)
    ConFiles = ploras_getconfiles(subjects,tasks(i),searchdir);
    roidata(:,i) = spm_summarise(ConFiles,struct('def','sphere', 'spec',region(4), 'xyz',region(1:3)'),stat);
end

    
    
    
    
    