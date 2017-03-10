# PLORAS MATLAB FUNCTIONS (for use at the FIL only)
---------------------------------------------------

This collection of MATLAB functions deals with all levels of processing of the PLORAS 3 fMRI paradigm. To be able to use these tools you have to add the ‘ploras_functions’ folder to your MATLAB search path. You can do this by typing the following commands in the MATLAB command window (you only have to do this once):

addpath Z:\PLORAS_fMRI_P3_01_SPM12_FEB2016\Scripts\ploras_functions;
savepath


## 1. Data import from CHARM

### ploras_data_import(src,dest,tweak) 

This function imports and extracts archives from specified source directories on 'Charm' into specified destination directories and converts images into '.nii' format. The function additionally sets up a folder structure in destination directory for functional and structural images and moves image files into their 
respective folders. 

Can be used in two ways:

1.	When processing only one subject: Two strings are required as input. The first string (SRC) is the source directory on 'Charm' and the second string (DEST) is the destination directory (both have to include the full path and have to be surrounded by inverted commas). 

Example: ploras_data_import('W:\trio\2013\Test1','Z:\PLORAS3_fMRI_01\Test1')

2.	When processing more than one subject: A single text file is required as input. The text file has to include two columns. The first column should contain the source directory names on 'Charm' and the second column should contain the corresponding destination directories (both have to include the full path). 

Example: 
ploras_data_import('folders.txt')

Requirement: 

Each destination directory must contain a run_info.txt file containing the run number - runtype correspondences separated by line breaks

(e.g. 
02 Fieldmap
03 Functional
04 Resting_state
05 Structural
...)

TWEAK (optional): Setting the third input argument to 1 allows to import data from a different scan as the "main" data of a subject (i.e. source and destination do not have to include the same scanID).

 
## 2. In-scanner behavioural performance

### ploras_update_behavioural_XL(subjects,searchdir,outlierthreshold)

This function updates the single-subject in-scanner behavioural performance Excel spreadsheets for the PLORAS 3 fMRI paradigm. This includes updating the accuracy and response times for the button-press tasks as well as updating the response times and durations for the speech production tasks. Additionally, all trial onsets are updated.

Required inputs:

SUBJECTS: either a single subject ID in inverted commas (e.g. 'PS0111’) or a text file containing a list of subject IDs (in a single column).

SEARCHDIR: The directory which contains the subject directories can be found (please provide the full path in inverted commas!).

OUTLIERTHRESHOLD (optional): specifies how many standard deviations from the mean constitutes 'speech' in the sound files. The default value is 2.

Requirement:
 
The script assumes that each subject folder has an existing behavioural Excel spreadsheet (at least a template).

### ploras_behavioural_summary_XL(subjects,tasks,searchdir,output,complex)

This script creates a summary Excel spreadsheet for the in-scanner behavioural performance for the PLORAS 3 fMRI paradigm. For each subject the output file contains a line with accuracy (i.e. percent correct 
trials) and mean response time (for correct trials) for each task.

Required inputs:

SUBJECTS: either a single subject ID in inverted commas (e.g. 'PS0111’) or a text file containing a list of subject IDs (in a single column).

TASKS: a numerical vector indicating the indices of the tasks of the P3 fMRI paradigm to be included. Examples are: 1:14 (for all tasks), [1 8 5 6] (for a subset of the tasks in a particular order). 

HINT: If you don’t know the index/indices of your task/s of interest, type ploras_fmri_conditions (for more info see below) in the MATLAB command line. This will list all tasks and their corresponding indices.

SEARCHDIR: The directory in which the subject directories can be found
(please provide the full path in inverted commas!).

OUTPUT: Name of the output Excel spreadsheet, e.g. '30_patients_summary'.
If you don't provide a full path, the file will be stored in the current
working directory.

COMPLEX (optional): If set to 0 (default) the output will contain a summary for the simple model of the fMRI paradigm. If set to 1, a summary for the complex model (i.e. 4 conditions per task) will be produced.

 
## 3. fMRI data preprocessing

## ploras_fmri_preprocessing(subjects,searchdir,TA)

This function performs the automatic preprocessing pipeline for the 
PLORAS 3 fMRI paradigm. Processing steps include realignment & unwarping of the functional MRI images, co-registration between functional and structural MRI images, normalisation to standard MNI space, and smoothing. 

Required inputs: 

SUBJECTS: the IDs of the subjects to include. Either a single subject ID
(e.g. 'PS0111') or a text file containing a list of subject IDs (in a single column). 

SEARCHDIR: the directory which includes the subject directories (please provide the full path in inverted commas!).

TA (optional): Set to 1, this will additionally run a time series analysis before and after the other preprocessing steps. This analysis visualises the mean and the variance of signal intensities per scan. The output of the time series analysis is printed to a .ps file in the subject-specific preprocessed_data folder.

Advanced options:

To switch off smoothing, make changes to the image resolution or the size of the smoothing kernel it is necessary to make manual edits to this function. To do this type ‘edit ploras_fmri_preprocessing.m’ in the Matlab command window. In the opened script go the ‘MANUAL EDITS’ section make the desired changes.

Before running the edited script make sure you save it under a different name (otherwise you will overwrite the default settings)!

 
## 4. Run fMRI first level analysis

### ploras_fmri_firstlevel(subjects,tasks,searchdir,complex)

This scripts performs the first level (i.e. subject-specific) analysis
for the PLORAS 3 fMRI paradigm. The resulting SPM.mat is stored in a newly created 'first_level' folder in the subject directory.

Required inputs:

SUBJECTS: either a single subject ID in inverted commas (e.g. 'PS0111’) or a text file containing a list of subject IDs (in a single column). 

TASKS: a numerical vector indicating the indices of the tasks of the P3 fMRI paradigm to be included. Examples are: 1:14 for all tasks, 1:5 or 
[1 2 3 4 5] for the first five tasks, or [1:3 5:14] when task no. 4 is missing.

HINT: If you don’t know the index/indices of your task/s of interest, type ploras_fmri_conditions (for more info see below) in the MATLAB command line. This will list all tasks and their corresponding indices.

SEARCHDIR: The directory in which the subject directories can be found
(please provide the full path in inverted commas!).

COMPLEX (optional): the analysis is either based on the simple model (0 = default) or the complex model (1).

Simple model:
•	4 regressors per task: instruction, correct, incorrect, missing
•	1 contrast per task: correct > fixation baseline
•	Results are stored in ‘\first_level\simple_model’

	Complex model:
•	7 regressors per task: instruction, 4 x correct, incorrect, missing
•	5 contrasts per task: all correct > fixation baseline + each condition > fixation baseline
•	Results are stored in ‘\first_level\complex_model’

 
## 5. Region of interest parameter estimate extraction

### ploras_roi_extraction(subjects,tasks,searchdir,region,stat)

This script extracts parameter estimates for specified task versus fixation baseline contrasts from voxels within a specified region of interest and calculates summary statistic across voxels for each specified subject. The resulting variable ROIDATA is a matrix of parameter estimates where each subject corresponds to a row and each task corresponds to a column. 

Required inputs:

SUBJECTS: Either a single subject ID in inverted commas (e.g. 'PS0111’) or a text file containing a list of subject IDs (in a single column). 

TASKS: a numerical vector indicating the indices of the tasks of the P3 fMRI paradigm to be included. Examples are: 1:14 for all tasks, 1:5 or [1 2 3 4 5] for the first 5 tasks, or [1:3 5:14] when task no. 4 is missing. 

Hint: If you don’t know the index/indices of your task/s of interest, type ploras_fmri_conditions (for more info see below) in the MATLAB command window. This will list all tasks and their corresponding indices.

SEARCHDIR: The directory which contains the subject directories can be found (please provide the full path in inverted commas!).

REGION: Region of interest. Can be either the name of an existing region of interest file (in inverted commas) or a numerical vector of length 4 specifying a spherical volume (i.e. [x, y, z, radius]).

STAT (optional): Summary statistic calculated across voxels included in the region of interest. Default is the first Eigenvariate of a single value decomposition algorithm (as is standard in SPM). An alternative is @mean which calculates the mean across voxels (which is not recommended however).

 
## 6. Motion parameters

### ploras_motion_parameters(subjects,tasks,searchdir,output)

This script extracts task-specific head motion parameters (path lengths) for the PLORAS 3 fMRI paradigm and saves them to an Excel spreadsheet.

Required inputs:

SUBJECTS: either a single subject ID in inverted commas (e.g. 'PS0111') or a text file containing a list of subject IDs (in a single column). If the same subject had more than one scan, the script will automatically detect this. This means that the ID does not have to be included more than once!

TASKS: a numerical vector indicating the indices of the tasks of the P3 fMRI paradigm to be included. Examples are: 1:14 (for all tasks), [1 8 5 6] (for a subset of the tasks in a particular order).

(HINT: If you don’t know the index/indices of your task/s of interest, type ploras_fmri_conditions in the MATLAB command line. This will list all tasks and their corresponding indices.)

SEARCHDIR: The directory in which the subject directories can be found
(please provide the full path in inverted commas!).

OUTPUT: Name of the output Excel spreadsheet, e.g. '30_patients_summary'.
If you don't provide a full path, the file will be stored in the current
working directory.

 
## 7. Utilities

### ploras_getconfiles(subjects,tasks,searchdir,excel,complex)

This function assembles the names (+ full path) of the contrast files for each subject for the indicated tasks, saves them in a output variable, and optionally saves them to an Excel spreadsheet in the current working directory 

Required inputs:

SUBJECTS: either a single subject ID in inverted commas (e.g. 'PS0111’) or a text file containing a list of subject IDs (in a single column). 

TASKS: a numerical vector indicating the indices of the tasks of the P3 fMRI paradigm to be included. Examples are: 1:14 (for all tasks), [1 8 5 6] (for a subset of the tasks in a particular order). If COMPLEX (see below) is set to 1, the function will save the contrast files for all 4 conditions of the tasks instead.

SEARCHDIR: The directory in which the subject directories can be found
(please provide the full path in inverted commas!).

EXCEL (optional): Name for output Excel file (without file extension) in 
inverted commas.

COMPLEX (optional): If set to 0 (default) the output will contain the contrast file names of the simple model analysis of the paradigm (i.e. 1 contrast per task: correct trials > fixation baseline). If set to 1, a summary for the complex model (i.e. 4 conditions per task) will be produced.

Examples:

ploras_getconfiles(‘20subjectx.txt’,1:5,’Z:\data\’,'Confiles_controls’,1)
confiles = ploras_getconfiles(‘20subjectx.txt’,1:5,’Z:\data\’) 

### ploras_fmri_conditions(index,complex)

This function can be used to get the names of the tasks and conditions
of the PLORAS 3 fMRI paradigm. 

When used without input arguments (i.e. ‘ploras_fmri_conditions’), the function will display all task names at the Matlab command window. 

INDEX (optional): a numerical vector indicating the indices of the tasks 
to be displayed (e.g. 1 or [1,5] or 1:5).

You can also save the task names by providing a output variable:

tasks = ploras_fmri_conditions

COMPLEX (optional): if this second input argument is set to 1, the function will additionally save the condition names of the (selected) tasks. In this case you have to specify 2 output variables:

[tasks conditions] = ploras_fmri_conditions(1:5,1)


