function ploras_update_behavioural_XL(subjects,searchdir,OutlierThreshold)
% PLORAS_UPDATE_BEHAVIOURAL_XL(SUBJECTS,SEARCHDIR,OUTLIERTHRESHOLD) updates 
% the single-subject in-scanner behavioural performance Excel spreadsheets for 
% the PLORAS 3 fMRI paradigm. This includes updating the accuracy and response 
% times for the button-press tasks as well as updating the response times and 
% durations for the speech production tasks. Additionally, all trial onsets 
% are updated.
%
% Required inputs:
%
% SUBJECTS: the IDs of the subjects to include. Either a single subject ID
% (e.g. 'PS0111') or a text file containing a list of subject IDs (in a single 
% column).
%
% SEARCHDIR: The directory which contains the subject directories can be found 
% (please provide the full path in inverted commas!).
%
% OUTLIERTHRESHOLD (optional): specifies how many standard deviations from 
% the mean constitutes 'speech' in the sound files. The default value is 3.
%
% Requirement:
% The script assumes that each subject folder has an existing behavioural Excel 
% spreadsheet (at least a template). 
%
% (c) Tom Hope / slightly adapted by Philipp L (last modified 01/12/2016)

warning off
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
% outlierthreshold
if nargin < 3
    OutlierThreshold = 3;
end

BufferS = 0.2; %RTs less than 200ms are anticipations.
BufferC = 0.1;

%% 1. Get all of the subjects and subject data file names
if(ischar(searchdir))
    AudioFolders = FindFolders(searchdir,'Audio');
else
    AudioFolders = searchdir;
end
SubjectFolders = cell(length(AudioFolders),1);
SubjectNames = cell(length(SubjectFolders),1);
SubjectXLs = cell(length(SubjectFolders),1);
for i=1:length(AudioFolders)
    SubjectNames(i) = {GetNameFromPath(AudioFolders{i})};
    SubjectFolders(i) = {substr(AudioFolders{i},0,(strfind(AudioFolders{i},'\behavioural\audio')-1))};
    BackSlash = strfind(SubjectFolders{i},'\');
    BackSlash = BackSlash(length(BackSlash));
    SubjectXLs(i) = {[SubjectFolders{i} '\behavioural\' substr(SubjectFolders{i},BackSlash,Inf), '.xlsx']};
end
BadFiles = {};

%% 2. For each subject, open the file, and work out what information needs
% to be added
for pf=1:length(SubjectXLs)
    [GoodData,BadWavs,BadRes] = GetWorkingFiles(AudioFolders{pf});
    if(~isempty(BadWavs))
        dlmcell(['BadWavs_' num2str(pf) '.txt'],BadWavs);
        dlmcell(['BadRes_' num2str(pf) '.txt'],BadRes);
    end    
    if(size(GoodData,1)>0)
        WavNames = GoodData(:,1);
        ResNames = GoodData(:,2);
        disp(['Subject ' num2str(pf) ': Found ' num2str(length(WavNames)) ' sound files.']);
        for j=1:length(WavNames)
            disp(['Processing file ' num2str(j) ' of folder ' num2str(pf) '...']);
            [Stream,Fs,~,~] = wavread(WavNames{j});
            High = 0.5;
            Low = 0.05;
            [hb,ha] = butter(10, High, 'low');
            [lb,la] = butter(10, Low, 'high');
            Stream = filtfilt(hb, ha, Stream);
            Stream = filtfilt(lb, la, Stream);
            
            Resolution = Fs / 1000;
            if(size(Stream,2) > 1)
                Maxes = mean(abs(Stream),1);
                MaxInd = find(Maxes == max(Maxes),1);
                Stream = Stream(:,MaxInd);
            end
            EventStartTimes = GetEventStartTimes(ResNames{1});
            if size(EventStartTimes,1) == 44
                EventStartTimes([1 12 23 34]) = [];
                EventStartTimes(:,2) = EventStartTimes(:,1) + 2500;
                Windows = floor(EventStartTimes .* Resolution);
                BaselineStream = Stream(Windows(20,2):Windows(21,1));
            elseif size(EventStartTimes,1) == 24
                EventStartTimes([1 7 13 19]) = [];
                EventStartTimes(:,2) = EventStartTimes(:,1) + 5000;
                Windows = floor(EventStartTimes .* Resolution);
                BaselineStream = Stream(Windows(10,2):Windows(11,1));
            end
            BaselineStream = BaselineStream((Fs*3):length(BaselineStream));
            BestWidth = FindWindow(BaselineStream,100,5000,100);
            %BestWidth = 3100;
            WindowLength = BestWidth;

            
            [Stream,Fs,~,~] = wavread(WavNames{j});
            High = 0.5;
            Low = 0.05;
            [hb,ha] = butter(10, High, 'low');
            [lb,la] = butter(10, Low, 'high');
             Stream = filtfilt(hb, ha, Stream);
             Stream = filtfilt(lb, la, Stream);
            Resolution = Fs / 1000;
            if(size(Stream,2) > 1)
                Maxes = mean(abs(Stream),1);
                MaxInd = find(Maxes == max(Maxes));
                Stream = Stream(:,MaxInd);
            end
            EventStartTimes = GetEventStartTimes(ResNames{j});
            EvenStartTimesWithInstructions = EventStartTimes - EventStartTimes(1);
            if size(EventStartTimes,1) == 44
                EventStartTimes([1 12 23 34]) = [];
                %ResultStruct(1) = EventStartTimes;
                EventStartTimes(:,2) = EventStartTimes(:,1) + 2500;
                Windows = floor(EventStartTimes .* Resolution);
                BaselineStream = Stream(Windows(10,2):Windows(11,1));
            elseif size(EventStartTimes,1) == 24
                EventStartTimes([1 7 13 19]) = [];
                %ResultStruct(1) = EventStartTimes;
                EventStartTimes(:,2) = EventStartTimes(:,1) + 5000;
                Windows = floor(EventStartTimes .* Resolution);
                BaselineStream = Stream(Windows(5,2):Windows(6,1));
            end
            
            SubjIDbegin = strfind(WavNames{j},'Pilot_') + 6;
            SubjIDbegin = SubjIDbegin(length(SubjIDbegin));
            SubjIDend = strfind(WavNames{j},'_');
            SubjIDend = SubjIDend(SubjIDend > SubjIDbegin);
            SubjIDend = SubjIDend(1);
            SubjID = substr(WavNames{j},SubjIDbegin-1,(SubjIDend-SubjIDbegin));

            TaskIDbegin = SubjIDend + 1;
            TaskIDend = strfind(WavNames{j},'_');
            TaskIDend = TaskIDend(TaskIDend > TaskIDbegin);
            TaskIDend = TaskIDend(length(TaskIDend)-1);
            TaskID = substr(WavNames{j},TaskIDbegin-1,(TaskIDend-TaskIDbegin));

            RTs = {};
            Durations = {};
            SpeakingOverTrialBoundary = false;
            for k=1:size(Windows,1)

                TickToTickSeries = Stream(Windows(k,1):Windows(k,2));
                %TickToTickSeries = TickToTickSeries(ceil(BufferS*Fs):length(TickToTickSeries));
                Means = WindowFilter(TickToTickSeries,WindowLength);
                [NewRT,NewDur,SpeakingOverTrialBoundary] = GetRTAndDur(Means,OutlierThreshold,SpeakingOverTrialBoundary,Resolution,BufferC,BufferS);
                
                RTs = [RTs ; {NewRT}];
                Durations = [Durations ; {NewDur}];
            end

            [pathname,RTfilesName] = fileparts(WavNames{j});
            RTfilesName = [pathname '\' RTfilesName '_RTs(ms).txt'];
            dlmcell(RTfilesName,[RTs Durations]);
            
            KeySheetName = TaskID;
            Num = GetTaskNum(TaskID);
            % Make the data into the right format
            if(length(RTs)==20)
                RTs = [{''} ; RTs(1:5) ; {''} ; RTs(6:10) ; {''} ; RTs(11:15) ; {''} ; RTs(16:20)];
                Durations = [{''} ; Durations(1:5) ; {''} ; Durations(6:10) ; {''} ; Durations(11:15) ; {''} ; Durations(16:20)];
                OnsetRange = 'C2:C25';
                RTRange = 'F2:F25';
                DurRange = 'J1:J25';
                StimDur = repmat({'5000'},5,1);
                StimDur = [{''} ; StimDur ; {''} ; StimDur ; {''} ; StimDur ; {''} ; StimDur];
                StimDurRange = 'D2:D25';
                Blanks = repmat({''},25,3);
                BlanksRange = 'K1:M25';
            else
                RTs = [{''} ; RTs(1:10) ; {''} ; RTs(11:20) ; {''} ; RTs(21:30) ; {''} ; RTs(31:40)];
                Durations = [{''} ; Durations(1:10) ; {''} ; Durations(11:20) ; {''} ; Durations(21:30) ; {''} ; Durations(31:40)];
                OnsetRange = 'C2:C45';
                RTRange = 'F2:F45';
                DurRange = 'J1:J45';
                StimDur = repmat({'2500'},10,1);
                StimDur = [{''} ; StimDur ; {''} ; StimDur ; {''} ; StimDur ; {''} ; StimDur];
                StimDurRange = 'D2:D45';
                Blanks = repmat({''},45,3);
                BlanksRange = 'K1:M45';
            end
            
           EvenStartTimesWithInstructions = num2cell(EvenStartTimesWithInstructions);
           WriteXLS(SubjectXLs{pf},Blanks,KeySheetName,BlanksRange);
           WriteXLS(SubjectXLs{pf},EvenStartTimesWithInstructions,KeySheetName,OnsetRange);
           WriteXLS(SubjectXLs{pf},RTs,KeySheetName,RTRange);
           WriteXLS(SubjectXLs{pf},[{'Resp. Duration'} ; Durations],KeySheetName,DurRange);
            if(~isempty(strfind(TaskID,'Pic')) || ~isempty(strfind(TaskID,'Read')) || ~isempty(strfind(TaskID,'Name')))
                if(ismember(Num,(1:4)))
                    StimDur = repmat({'2500'},5,1);
                    StimDur = [{''} ; StimDur ; {''} ; StimDur ; {''} ; StimDur ; {''} ; StimDur];
                    StimDurRange = 'D2:D25';
                elseif(ismember(Num,[6 8 9 11]))
                    StimDur = repmat({'1500'},10,1);
                    StimDur = [{''} ; StimDur ; {''} ; StimDur ; {''} ; StimDur ; {''} ; StimDur];
                    StimDurRange = 'D2:D45';
                end
                WriteXLS(SubjectXLs{pf},StimDur,KeySheetName,StimDurRange);
            end
        end
    end
    % Now we have to get the assessment trial data
    MatchingTaskFolder = [SubjectFolders{pf} '\behavioural\matching_tasks'];
    MatchingResFiles = dir([MatchingTaskFolder '\*.txt']);
    MatchingResFiles = {MatchingResFiles.name}';
    for i=1:length(MatchingResFiles)
        TaskID = GetTaskIDFromResFile(MatchingResFiles{i});
        u = strfind(TaskID,'_');
        u = u(length(u)-1);
        KeySheetName = substr(TaskID,0,u-1);
        data = textread([MatchingTaskFolder '\' MatchingResFiles{i}],'%s','delimiter','\t');
        % Should get a 274,1 cell array (if the file is full)
        GoodFile = 1;
        if(size(data,1)==274)
            EventOnsetInds = (21:11:274);
            AccInds = (16:11:274);
            RTInds = (17:11:274);
            Onsets = data(EventOnsetInds);
            Acc = data(AccInds);
            RT = data(RTInds);
        elseif(size(data,1)==201)
            EventOnsetInds = (16:8:201);
            AccInds = (14:8:201);
            RTInds = (15:8:201);
            Onsets = data(EventOnsetInds);
            Init = str2num(Onsets{1});
            for a=1:length(Onsets)
                Onsets(a) = {num2str(str2num(Onsets{a})-Init)};
            end
            
            Acc = data(AccInds);
            RT = data(RTInds);
        else
            GoodFile = 0;
        end
        if(GoodFile > 0)
            for a=1:length(Acc)
                if(strcmp(Acc{a},'-1'))
                    Acc(a)={''};
                else
                    Acc(a)={num2str(str2num(Acc{a})+1)};
                end
                if(strcmp(RT{a},'-1')),RT(a)={''};end
                if(~isempty(strfind(RT{a},'m'))||~isempty(strfind(RT{a},'n')))
                    RT(a)={''};
                    Acc(a) = {'0'};
                end
            end
            if(length(Onsets)<24) % incomplete file
                Difference = 24-length(Onsets);
                FinalStr = num2str(25-Difference);
                OnsetRange = ['C2:C' FinalStr];
                RTRange = ['F2:F' FinalStr];
                AccRange = ['E2:E' FinalStr];
            else
                OnsetRange = 'C2:C25';
                RTRange = 'F2:F25';
                AccRange = 'E2:E25';
            end

            StimDur = repmat({'2500'},5,1);
            StimDur = [{''} ; StimDur ; {''} ; StimDur ; {''} ; StimDur ; {''} ; StimDur];
            StimDurRange = 'D2:D25';
            Blanks = repmat({''},25,3);
            BlanksRange = 'K1:M25';
            WriteXLS(SubjectXLs{pf},Blanks,KeySheetName,BlanksRange);
            WriteXLS(SubjectXLs{pf},Acc,KeySheetName,AccRange);
            WriteXLS(SubjectXLs{pf},RT,KeySheetName,RTRange);
            WriteXLS(SubjectXLs{pf},Onsets,KeySheetName,OnsetRange);
            if(~isempty(strfind(TaskID,'Pic')))
                WriteXLS(SubjectXLs{pf},StimDur,KeySheetName,StimDurRange);
            end
        else
            Blanks = repmat({''},25,3);
            BlanksRange = 'K1:M25';
            WriteXLS(SubjectXLs{pf},Blanks,KeySheetName,BlanksRange);
            WriteXLS(SubjectXLs{pf},{'Bad Results File'},KeySheetName,'K1:K1');
            dlmcell(['BadRes_' num2str(pf) '.txt'],[SubjectXLs{pf} KeySheetName]);
        end
    end
        
end
end

function WriteXLS(XLS,Data,SheetName,Range)
    Attempts = 0;
    Done = false;
    while(Attempts < 5 && Done == false)
        try
            xlswrite(XLS,Data,SheetName,Range);
            Done = true;
        catch ME
            disp(ME.message)
            disp(XLS)
            disp(Data)
            disp([SheetName ': ' Range])
        end
        Attempts = Attempts + 1;
    end
end

function Num = GetTaskNum(TaskID)
    u = strfind(TaskID,'_');
    Num = substr(TaskID,0,u(1)-1);
    Num = str2num(Num);
end

function AudioFolderNames = FindFolders(RootFolder,Mode)
    % First, we want to find all of the 'audio' folders
    AudioFolders = fdir([RootFolder '\**\behavioural\audio\']);
    AudioFolderNames = {};
    for i=1:length(AudioFolders.folders)
        if(~isempty(strfind(AudioFolders.folders{i},'audio')))
            AudioFolderNames = [AudioFolderNames ; AudioFolders.folders(i)];
        end
    end
    if(strcmp(Mode,'Subjects'))
        for i=1:length(AudioFolderNames)
            AudioFolderNames(i) = {substr(AudioFolderNames{i},0,(strfind(AudioFolderNames{i},'\behavioural\audio')-1))};
        end
    end
    disp(['Found ' num2str(length(AudioFolderNames)) ' folders.']);
end

function SubjName = GetNameFromPath(folder)
    backslashes = strfind(folder,'\');
    backslashes = backslashes(backslashes<strfind(folder,'\behavioural\audio'));
    backslashes = backslashes(length(backslashes));
    keybit=substr(folder,backslashes,strfind(folder,'\behavioural\audio')-(backslashes+1));
    u=strfind(keybit,'_');
    startu=u(5);
    endu=u(6);
    SubjName=substr(keybit,startu,endu-(startu+1));
end

function [GoodData,BadWavs,BadRes] = GetWorkingFiles(AudioFolder)

    [GoodWavs,BadWavs] = GetGoodWavs(AudioFolder);
    [GoodRes,BadRes] = GetGoodRes(AudioFolder);
    
    GoodData = {}; GoodInd = 0;
    for i=1:length(GoodWavs)
        ThisWav = GoodWavs{i};
        
        SubjIDbegin = strfind(ThisWav,'Pilot_') + 6;
        SubjIDbegin = SubjIDbegin(length(SubjIDbegin));
        SubjIDend = strfind(ThisWav,'_');
        SubjIDend = SubjIDend(SubjIDend > SubjIDbegin);
        SubjIDend = SubjIDend(1);
        
        TaskIDbegin = SubjIDend;
        TaskIDend = strfind(ThisWav,'_');
        TaskIDend = TaskIDend(TaskIDend > TaskIDbegin);
        TaskIDend = TaskIDend(3);
        TaskID = substr(ThisWav, TaskIDbegin, (TaskIDend-(TaskIDbegin+1)));

        j=1;
        while(j<=length(GoodRes) && isempty(strfind(GoodRes{j},TaskID)))
            j=j+1;
        end
        if(j<=length(GoodRes))
            GoodInd = GoodInd + 1;
            GoodData(GoodInd,1) = {ThisWav};
            GoodData(GoodInd,2) = GoodRes(j);
        end
    end
end

function [GoodWavs,BadWavs] = GetGoodWavs(AudioFolder)
    WavFiles = rdir([AudioFolder '\*.wav']);
    WavNames = {WavFiles.name}';
    disp(['Found ' num2str(length(WavNames)) ' sound files.']);
    % Exclude all of the 'practice, 'scanner practice', and 'test' files
    % from both lists
    Excluder = ones(length(WavNames),1);
    for j=1:length(WavNames)
        if(~isempty(strfind(WavNames{j},'Scanner Practice')) || ...
                ~isempty(strfind(WavNames{j},'Practice')) || ...
                ~isempty(strfind(WavNames{j},'5_Aud_Sem_Ass_2_Ob')) || ...
                ~isempty(strfind(WavNames{j},'1_Pic_Sem_Ass_2_Ob')))
            Excluder(j) = 0;
        end
    end
    WavNames = WavNames(logical(Excluder));
    WavNames = sort(WavNames);
    Selector = true(length(WavNames),1);
    for i=1:length(WavNames)
        try
            temp = wavread(WavNames{i});
        catch
            Selector(i) = false;
        end
    end
    BadWavs = WavNames(~Selector);
    GoodWavs = WavNames(Selector);
end

function [GoodRes,BadRes] = GetGoodRes(AudioFolder)
    ResFolder = strrep(AudioFolder,'audio','speaking_tasks');
    ResFiles = rdir([ResFolder '\*.txt']);
    ResNames = {ResFiles.name}';
    disp(['Found ' num2str(length(ResNames)) ' sound files.']);
    % Exclude all of the 'practice, 'scanner practice', and 'test' files
    % from both lists
    Excluder = ones(length(ResNames),1);
    for j=1:length(ResNames)
        if(~isempty(strfind(ResNames{j},'Scanner Practice')) || ...
                ~isempty(strfind(ResNames{j},'Practice')) || ...
                ~isempty(strfind(ResNames{j},'5_Aud_Sem_Ass_2_Ob')) || ...
                ~isempty(strfind(ResNames{j},'1_Pic_Sem_Ass_2_Ob')))
            Excluder(j) = 0;
        end
    end
    ResNames = ResNames(logical(Excluder));
    ResNames = sort(ResNames);
    Selector = true(length(ResNames),1);
    for i=1:length(ResNames)
        EventStartTimes = GetEventStartTimes(ResNames{i});
        if(isempty(EventStartTimes))
            Selector(i) = false;
        end
    end
    BadRes = ResNames(~Selector);
    GoodRes = ResNames(Selector);
end

function TaskID = GetTaskIDFromResFile(FileName)
    u = strfind(FileName,'_');
    Start = u(8);
    End = u(length(u)-1) - 1;
    Len = End - Start;
    TaskID = substr(FileName,Start,Len);
end

function EventStartTimes = GetEventStartTimes(ResFile)
    
    EventStartTimes = []; 
    try
        temp = importdata(ResFile);
    catch
        return; % File not found
    end
    
    try
        EventStartTimes = temp.data(:,3);
    catch
        try
            fid = fopen(ResFile);
            FileData = textscan(fid,'%s %s %s %d %d %[^\n]','headerlines',1);
            EventStartTimes = cell2mat(FileData(1,5));
        catch
            return;
        end
    end
end

function BestWidth = FindWindow(Series,Initial,Max,Increment)
BestWidth = 0; BestDev = Inf; GraphDev = [];
FirstMinFound = 0; Width = Initial;
while(Width <= Max && FirstMinFound == 0)
    End = 0; Start = 1;
    Means = ones(length(Series),1).*-1;
    while(End < length(Series))
        End = Start + Width;
        Means(Start,1) = mean(abs(Series(Start:End)));
        Start = Start + 1;
        %disp(['W=' num2str(Width) '; S=' num2str(Start)]);
    end
    Means(Means<0) = [];
    MaxDeviation = max(abs(Means - mean(Means)));
    if(MaxDeviation < BestDev)
        BestDev = MaxDeviation;
        BestWidth = Width;
    elseif(MaxDeviation > BestDev)
        FirstMinFound = 1;
    end
    GraphDev = [GraphDev ; MaxDeviation];
    plot(GraphDev); drawnow; pause(1);
    Width=Width+Increment;
end
end

function Filtered = WindowFilter(Series,WindowLength)
Start = 1; Means = [];
while((Start + (WindowLength-1)) <= length(Series))
    End = Start + (WindowLength-1);
    Means(Start,1) = mean(abs(Series(Start:End)));
    Start = Start + 1;
end
pause(.001);
Filtered = zscore(Means);
end

function [NewRT,NewDur,SpeakingOverTrialBoundary] = GetRTAndDur(Means,OutlierThreshold,SpeakingOverTrialBoundary,Resolution,BufferC,BufferS)
try
AboveThreshold = find(Means>=OutlierThreshold);
AtBaseline = find(Means<0);
BeforeFirst = AtBaseline(AtBaseline<AboveThreshold(1));
FirstOnset = BeforeFirst(length(BeforeFirst));
if(SpeakingOverTrialBoundary)
    NewRT = '';
    NewDur = '';
    SpeakingOverTrialBoundary = false;
    LastOnset = AboveThreshold(length(AboveThreshold));
    if(LastOnset >= length(Means))
        SpeakingOverTrialBoundary = true;
    end
else
    NewRT = (FirstOnset / Resolution) + (BufferS * 1000);
    LastOnset = AboveThreshold(length(AboveThreshold));
    if(LastOnset < length(Means))
        ReturnPostSpeech = AtBaseline(AtBaseline > LastOnset);
        ReturnPostSpeech = ReturnPostSpeech(1);
        NewDur = (ReturnPostSpeech - FirstOnset) / Resolution;
        SpeakingOverTrialBoundary = false;
    else
        NewDur = {''};
        SpeakingOverTrialBoundary = true;
    end
end
catch
    NewRT = '';
    NewDur = '';
    SpeakingOverTrialBoundary = false;
end
end
