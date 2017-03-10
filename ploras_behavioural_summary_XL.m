function ploras_behavioural_summary_XL(subjects,tasks,searchdir,output,complex)
% PLORAS_BEHAVIOURAL_SUMMARY_XL(SUBJECTS,TASKS,SEARCHDIR,OUTPUT,COMPLEX) 
% creates a summary Excel spreadsheet for the in-scanner behavioural performance 
% for the PLORAS 3 fMRI paradigm. For each subject the file contains a line 
% with accuracy (i.e. percent correct trials) and mean response time (for 
% correct trials) for each task.
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
% COMPLEX (optional): If set to 0 (default) the output will contain a summary for
% the simple model of the paradim. If set to 1, a summary for the complex
% model (i.e. 4 conditions per task) will be produced.
%
% Philipp Ludersdorfer (last modified 29/11/2016)


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

% task list
if(isnumeric(tasks))
    TaskList = ploras_fmri_conditions(tasks);
end
% output
if ~ischar(output)
    error('Invalid output file name!')
end
% simple version (skip data for 4 conditions within each task)
if nargin < 5
    complex = 0;
end

SummaryArray = {}; wrngs = {}; wrn = 1;

%% LOOP OVER SUBJECTS/FILES
for i=1:length(FinalFolders)
    % Open behavioural XL spreadsheet of subject
    disp('# # # # # # # # # # # # # # # # # # # # # # # # # #')
    disp(['Working on subject ' FinalSubjects{i} ' ' num2str(i) '/' num2str(length(FinalFolders))])
    disp('# # # # # # # # # # # # # # # # # # # # # # # # # #')
    xlsdir = dir(fullfile(searchdir,FinalFolders{i},'behavioural','*.xlsx'));
    try
        KeyXL = fullfile(searchdir,FinalFolders{i},'behavioural',xlsdir(1).name);
        CurrentXL = importdata(KeyXL);
    catch
        wrngs{wrn} = [FinalSubjects{i} ': No existing behavioural excel file'];
        wrn = wrn + 1;
        continue
    end
    % LOOP OVER TASKS/EXCEL SHEETS
    for t=1:length(TaskList)
        disp(['Task ' num2str(t) '/' num2str(length(TaskList))]);
        KeySheet = ['x' TaskList{t}];
        Text = CurrentXL.textdata.(KeySheet);
        Data = CurrentXL.data.(KeySheet);
        RTs = Data(:,4);
        Acc = Data(:,3);
        Code = Data(:,6);
        Code2 = Data(:,7);
        RTs = RTs(~isnan(Acc));
        Code = Code(~isnan(Acc));
        Code2 = Code2(~isnan(Acc));
        Codes = [Code,Code2];
        uCodes = unique(Codes,'rows');
        Acc = Acc(~isnan(Acc));
        Underscore = strfind(KeySheet,'_');
        Underscore = Underscore(1);
        Num = str2num(substr(KeySheet,1,Underscore-2)); %#ok<*ST2NM>
        if(isempty(Num))
            warn('Num not working');
            disp(Num);
        elseif(Num == 2 || Num == 4)
            CorrectOnly = Acc >= 4;
        else
            CorrectOnly = Acc == 2;
        end
        % Get summary stats (simple model)
        [minval_All,maxval_All,meanval_All,medianval_All,accrate_All] = GetSummaryStats(RTs,CorrectOnly);
        WriteData = [accrate_All,meanval_All,medianval_All,minval_All,maxval_All];
        % Get summary stats (complex model)
        if complex
            [minval_COD,maxval_COD,meanval_COD,medianval_COD,accrate_COD] = GetSummaryStats(RTs,CorrectOnly,Code,Code2,uCodes);
            for w=1:length(accrate_COD)
                WriteData = [WriteData,accrate_COD(w),meanval_COD(w),medianval_COD(w),minval_COD(w),maxval_COD(w)];
            end
        end
        % Setup SummaryArray
        if(isempty(SummaryArray))
            SummaryArray = SetupSummaryArray(FinalSubjects,TaskList,uCodes,complex);
        end
        % Format stats for writing into Excel File
        WriteRow = 2+i;
        WriteCols = 6+t;
        for c=2:length(WriteData)
            WriteCols = [WriteCols (WriteCols(length(WriteCols))+length(TaskList))];
        end
        for c=1:length(WriteCols)
            SummaryArray(WriteRow,WriteCols(c)) = {WriteData(c)};
        end
        [Native,DomHand,Lesion,RespSide,Scanner] = GetMetaData(FinalFolders{i});
        SummaryArray(WriteRow,2:6) = [{Native} {DomHand} {Lesion} {RespSide} {Scanner}];
    end
    disp('Done.')
end
% Write everything to Excel spreadsheet
xlswrite([output '.xlsx'],SummaryArray);

%% DONE
disp('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #')
disp('The behavioural summary spreadheet has been successfully created!')
disp('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #')
disp(' ')
if wrn~= 1
    disp('The following errors occurred:')
    for w = 1:length(wrngs)
        disp(wrngs{w})
    end
end

end

function SummaryArray = SetupSummaryArray(SubjList,TaskList,uCodes,complex)

    SummaryArray={};
    SummaryArray(2,1:6) = {'SUBJECT' 'NATIVE' 'DOMHAND' 'LESION' 'RESPHAND' 'SCANNER'};
    SummaryArray(3:(2+length(SubjList)),1) = SubjList;
    Start = 7; End = (Start-1) + length(TaskList);
    SummaryArray(1,Start:End) = repmat({'ACCURACY'},1,length(TaskList));
    SummaryArray(2,Start:End) = TaskList';
    Start = End + 1; End = (Start-1) + length(TaskList);
    SummaryArray(1,Start:End) = repmat({'MEAN RTs Correct Only'},1,length(TaskList));
    SummaryArray(2,Start:End) = TaskList';
    Start = End + 1; End = (Start-1) + length(TaskList);
    SummaryArray(1,Start:End) = repmat({'MEDIAN RTs Correct Only'},1,length(TaskList));
    SummaryArray(2,Start:End) = TaskList';
    Start = End + 1; End = (Start-1) + length(TaskList);
    SummaryArray(1,Start:End) = repmat({'MIN RTs Correct Only'},1,length(TaskList));
    SummaryArray(2,Start:End) = TaskList';
    Start = End + 1; End = (Start-1) + length(TaskList);
    SummaryArray(1,Start:End) = repmat({'MAX RTs Correct Only'},1,length(TaskList));
    SummaryArray(2,Start:End) = TaskList';
    
    if complex
        for i=1:length(uCodes)
            Start = End + 1; End = (Start-1) + length(TaskList);
            CodeLabel = ['Code_' num2str(uCodes(i,1)) '_' num2str(uCodes(i,2)) '_'];
            SummaryArray(1,Start:End) = repmat({[CodeLabel 'ACCURACY']},1,length(TaskList));
            SummaryArray(2,Start:End) = TaskList';
            Start = End + 1; End = (Start-1) + length(TaskList);
            SummaryArray(1,Start:End) = repmat({[CodeLabel 'MEAN RTs Correct Only']},1,length(TaskList));
            SummaryArray(2,Start:End) = TaskList';
            Start = End + 1; End = (Start-1) + length(TaskList);
            SummaryArray(1,Start:End) = repmat({[CodeLabel 'MEDIAN RTs Correct Only']},1,length(TaskList));
            SummaryArray(2,Start:End) = TaskList';
            Start = End + 1; End = (Start-1) + length(TaskList);
            SummaryArray(1,Start:End) = repmat({[CodeLabel 'MIN RTs Correct Only']},1,length(TaskList));
            SummaryArray(2,Start:End) = TaskList';
            Start = End + 1; End = (Start-1) + length(TaskList);
            SummaryArray(1,Start:End) = repmat({[CodeLabel 'MAX RTs Correct Only']},1,length(TaskList));
            SummaryArray(2,Start:End) = TaskList';
        end
    end
end

function [Native,DomHand,Lesion,RespSide,Scanner] = GetMetaData(folder)
u=strfind(folder,'_');
Native=substr(folder,u(1),u(2)-(u(1)+1));
DomHand=substr(folder,u(2),u(3)-(u(2)+1));
if strcmp(folder(1),'P')
    Lesion=substr(folder,u(3),u(4)-(u(3)+1));
    RespSide=substr(folder,u(4),u(5)-(u(4)+1));
    if length(u) == 11
        Scanner=substr(folder,u(10),u(11)-(u(10)+1));
    else
        Scanner=substr(folder,u(9),u(10)-(u(9)+1));
    end
else
    Lesion='';
    RespSide=substr(folder,u(3),u(4)-(u(3)+1));
    Scanner=substr(folder,u(8),u(9)-(u(8)+1));
end
end

function [minval,maxval,meanval,medianval,accrate] = GetSummaryStats(RTs,Acc,Coding,Coding2,uCodes)
[minval, maxval, meanval, medianval, accrate] = deal(-1);
if(nargin == 1)
    minval = min(RTs);
    maxval = max(RTs);
    meanval = mean(RTs);
    medianval = median(RTs);
    accrate = -1;
elseif(nargin==2)
    S = Acc & ~isnan(RTs);
    if(~isempty(find(S,1)))
        minval = min(RTs(S));
        maxval = max(RTs(S));
        meanval = mean(RTs(S));
        medianval = median(RTs(S));
    end
        accrate = (length(find(Acc)) / length(Acc)) * 100;
elseif(nargin == 5)
    for c=1:length(uCodes)
        S = Acc & ~isnan(RTs) & Coding==uCodes(c,1) & Coding2==uCodes(c,2);
        AccS = Acc & Coding==uCodes(c,1) & Coding2==uCodes(c,2);
        Sa = Coding==uCodes(c,1) & Coding2==uCodes(c,2);
        if(~isempty(find(S,1)))
            minval(c,1) = min(RTs(S));
            maxval(c,1) = max(RTs(S));
            meanval(c,1) = mean(RTs(S));
            medianval(c,1) = median(RTs(S));
            accrate(c,1) = (length(find(AccS)) / length(AccS)) * 100;
        else
            minval(c,1) = -1;
            maxval(c,1) = -1;
            meanval(c,1) = -1;
            medianval(c,1) = -1;
            accrate(c,1) = (length(find(AccS)) / length(AccS)) * 100;
        end
    end
end
end