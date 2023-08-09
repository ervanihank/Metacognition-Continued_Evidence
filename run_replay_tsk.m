% script to run Replay Task with drifting gabors stimuli

% each trial the observer is shown an arrage of gabors, whose phase is
% updated to appear as if drifting in a particular direction
% the observer must decide as quickly as possible if the gabors are
% drifting leftward or rightward. In one trial type, the stimuli are 
% displayed until the response is entered. In another, the stimuli are
% displayed for a fixed duration, which is determined based on the previous
% duration that the participant chose to view those stimuli for.

%--------------------------------------
%       Script Outline
%--------------------------------------

% 1. Set up
%   a. Directory
%   b. Participant number
%   c. Data storage
%   d. Parameters and Conditions
%   e. Screen, Sound, keyboard
%   f. Stimuli
% 2. Trials
%   b. Trial Loop
%       i. iti
%       ii. Stimulus
%       iii. response
%   c. data sorage
% 3. Analysis
% 4. Clean up

%-------------------------------------
%       Input
%-------------------------------------

% freeDataFileName:     the name of the data file from the free task

% stimParams are loaded from the saved file


%-----------------------------------
%       Output
%-----------------------------------

% p:             parameters
% data:          data 
% stimParams:    stimulus parameters for regenerating stimuli


%-----------------------------------
%       1. Set up
%------------------------------------


function data = run_replay_tsk(freeDataFolderName, freeDataFileName)

Screen('Preference', 'SkipSyncTests', 1);

try

commandwindow;

experimentName = 'contEvExp';
taskName = 'replayTsk';

% Directory

currentDir = cd;
experimentDir = [filesep experimentName];

% check that we are in the correct directory
correctDir = strcmp(currentDir(end-length(experimentDir)+1:end), experimentDir);

if correctDir
    addpath(genpath(currentDir))
else
    fprintf(['\n\n YOU ARE IN THE WRONG DIRECTORY!!! \n\n Please move to ' experimentDir '\n\n'])
    return
end


% load the previous data
data_folder = [cd filesep 'data_' experimentName];

% the folder name is the participant id, which is the first characters
% before the '_'
% ptchar = 1;
% charcomp = 0;
% while ~charcomp
%     if strcmp(freeDataFileName(ptchar),'_')
%         charcomp = 1;
%     else
%         ptchar = ptchar+1;
%     end
% end

prevdat = load([data_folder filesep freeDataFolderName filesep freeDataFileName]);

p = prevdat.p;
dataF = prevdat.data;

% load the stimulus parameters
load([cd filesep 'stimParamsFree.mat'])

% Data storage

dataSavePath = [data_folder filesep freeDataFolderName];

% previous data files
dataFiles = ls(dataSavePath);
if ~isempty(dataFiles)
    fprintf(['\n\nHere are the existing data files for ', p.participant, ':\n\n']);
    disp(dataFiles);
end

allFiles = dir(dataSavePath);
validFiles = 0;

for f = 1:length(allFiles)
    if ~isempty(strfind(allFiles(f).name, experimentName))
        validFiles = validFiles + 1;
        runFiles{validFiles, :} = allFiles(f).name;
    end
end

% session number

fprintf('\n');
for f = 1:validFiles
    fprintf(['(', int2str(f), ')  ', runFiles{f}, '\n']);
end
session = ( input('Please enter the session number:     ') );

sessiontype = ['S' num2str(session)];

fprintf('\n-------------------------------\n');
fprintf([' THIS RUN IS SESSION ' int2str(session) '\n']);
fprintf('-------------------------------\n\n');


% Get date and time for trial
time_stamp = datestr(now, 31);       % get current date and time as string
time_stamp = time_stamp(1:end-3);    % chop seconds off
time_stamp(end-[2 5]) = '-';         % replace space and colon with dashes

% Name files
dataSaveFileName = strcat(p.participant, '_', experimentName, '_', taskName, '_', sessiontype , '_', time_stamp);
fullSaveFileName = [dataSavePath filesep dataSaveFileName '.mat'];  % final filename with full path
disp(['New data filename is: ', fullSaveFileName]);    
    

% update parameters

p.session = session;

% set up random seed to be able to regenerate experiment
p.startTime = time_stamp;
p.dataFile = fullSaveFileName;
p.taskName = taskName;

p.initRandSeed = GetSecs;
rng(p.initRandSeed);

% conditions

% we will have: 
% 3 x free trials
% 3 x less trials
% 2 x more -
% 2 x more /
% 2 x more +
p.n_rep = 3+3+2+2+2;
p.n_trials = p.n_rep*p.n_indiv_trials;

% condition labels
condlabelsu = [1,1,1,2,2,2,3,3,4,4,5,5];

% whether the trial has a fixed duration
fixTime = [0,0,0,1,1,1,...
    1,1,1,1,1,1];

% condition labels

condLabels = sort(repmat(condlabelsu,1,p.n_indiv_trials));
fixTimes = sort(repmat(fixTime,1,p.n_indiv_trials));

trialConditions = zeros(p.n_trials,5);
trialConditions(:,1:3) = repmat(p.conditionMat,p.n_rep,1);
trialConditions(:,4) = condLabels;
trialConditions(:,5) = fixTimes;

% shuffle the order
shuffInds = randperm(p.n_trials);
p.trialConditions = trialConditions(shuffInds,:);


% calculate the durations for each trial
freeDec = zeros(p.n_indiv_trials,1);

for ti = 1:p.n_indiv_trials
    
    thesets = dataF.data(:,3) == ti;
    
    medrt = median(dataF.data(thesets,8));
    
    freeDec(ti) = medrt - 0.2;
    
end

freeDec(freeDec<=0.05) = 0.05;

lessDur = freeDec - 0.3;
moreDur = freeDec + 0.3;

% check for any too short
lessDur(lessDur<0.1) = 0.1;

% or waaaay too long
freeDec(freeDec>2) = 2;

% convert to frames
freeDecFlips = round(freeDec/p.flipInterval);
lessFlips = round(lessDur/p.flipInterval);
moreFlips = round(moreDur/p.flipInterval);


% generate the stimulus parameters for the more conditions
[stimParamsN, stimParamsS, stimParamsC] = genStimParamsReplay(stimParams, freeDecFlips, moreFlips-freeDecFlips, p.n_gabors, p.rangeMin, p.p_noi, p.win_size,p.lifetime,p.drift_speed);

% data
data.headers = {'block'; 'Trial'; 'stimID'; 'range'; 'trueDir'; 'ReplayCondition'; 'Response'; 'Correct'; 'RT'; 'Confidence'; 'ConfidenceRT'};
data.data = zeros(p.n_trials, length(data.headers));
data.prevResps = cell(p.n_trials,1);



% Screen, Sound, Keyboard
if ~p.dummyMode
    [window, rect] = PsychImaging('OpenWindow', p.screenNo, p.grey);
else
    [window, rect] = PsychImaging('OpenWindow', p.screenNo, p.grey,[0,800;0,800]);
end

Screen('BlendFunction',window, 'GL_ONE','GL_ZERO');

% Center of the screen
[p.center_x, p.center_y] = RectCenter(rect);


% check flip interval
p.flipInterval = Screen('GetFlipInterval', window);
refreshRate = 1 / p.flipInterval;
if round(refreshRate*10) ~= round(p.desiredRefreshRate*10)
    fprintf('\nReal refresh rate is: %2.2f Hz\n', refreshRate);
    p.screenRefreshRateActual = refreshRate;
end

topPriorityLevel = MaxPriority(window);

HideCursor;

Screen('FillRect', window, p.grey, rect); % clear screen
startvbl = Screen('Flip', window);
 
% keyboard set up

%keyboard set up
k.kbDev = -1;
KbName('UnifyKeyNames');
% only check the keys we need - left arrow, right arrow and escape
% plus confidence - 1-4
k.escapeKey = KbName('ESCAPE');
k.right = KbName('rightArrow');
k.left = KbName('leftArrow');
k.one = KbName('1!');
k.two = KbName('2@');
k.three = KbName('3#');
k.four = KbName('4$');
k.space = KbName('space');

ListenChar(-1);

% List for KbQueue
k.queueList = zeros(1,256);
k.queueList([k.escapeKey, k.right, k.left, k.space, k.one, k.two, k.three, k.four]) = 1;

KbQueueCreate(k.kbDev,k.queueList);

% stimuli

% create the procedural gabor
gabor_tex = CreateProceduralGabor(window, p.gabor_text_dim, p.gabor_text_dim,[],p.backgroundgabor_text , p.disable_norm, p.pre_contrast_multiplier);

%--------------------------------------
%          2. trials
%---------------------------------------


% loop through blocks where the participant is given a break

% block look
blocks = 1:(p.n_trials/p.nBlocks):p.n_trials;

for b = 1:p.nBlocks
    
    % reset queue
    KbReleaseWait(k.kbDev);
    KbQueueFlush(k.kbDev);
    
    % start new queue
    KbQueueStart(k.kbDev);
    
    thisBlock = blocks(1,b);
    
    blockMessage = ['This is block number ' num2str(b) ' of ' num2str(p.nBlocks) '\n\n' ' \n\n Press space to begin!'];
    
    DrawFormattedText(window,blockMessage,'center','center',p.black);
    
    Screen('Flip',window); %display on screen
    
    pressed = 0;
    doContinue = 0;
    while ~doContinue
        while ~pressed
            [pressed, keys] = KbQueueCheck(k.kbDev);
        end
        
        key = min(find(keys,1,'first'));
        if key == k.escapeKey
            fprintf('User Quit')
            clean_up
            return
        elseif key == k.space
            doContinue = 1;
        else
            pressed = 0;
        end
        
    end
    
    KbQueueStop(k.kbDev);
    
    for t = thisBlock:(thisBlock + (p.n_trials/p.nBlocks) -1)
        
        
        stimTrialConds = p.trialConditions(t,:);
        
        switch stimTrialConds(4)
            
            case 1 % free trial
                theseOrientations = stimParams(stimTrialConds(1)).orientation';
                thesePhases = stimParams(stimTrialConds(1)).phases';
                thisMaxFlips = p.maxFlips;
            case 2 % less trial
                thisMaxFlips = lessFlips(stimTrialConds(1));
                theseOrientations = stimParams(stimTrialConds(1)).orientation';
                thesePhases = stimParams(stimTrialConds(1)).phases(1:thisMaxFlips,:)';
                
            case 3 % more minus trial
                theseOrientations = stimParamsC(stimTrialConds(1)).orientation';
                thesePhases = stimParamsC(stimTrialConds(1)).phases';
                thisMaxFlips = moreFlips(stimTrialConds(1));
            case 4 % more same trial
                theseOrientations = stimParamsN(stimTrialConds(1)).orientation';
                thesePhases = stimParamsN(stimTrialConds(1)).phases';
                thisMaxFlips = moreFlips(stimTrialConds(1));
            case 5 % more plus trial
                theseOrientations = stimParamsS(stimTrialConds(1)).orientation';
                thesePhases = stimParamsS(stimTrialConds(1)).phases';
                thisMaxFlips = moreFlips(stimTrialConds(1));
        end
                
        
        
        % present the stimulus and get the response
        
        [response, rt, conf, confRT, prevResps] = presentTrialReplay(window, p, k, gabor_tex, theseOrientations, thesePhases, stimTrialConds(5),thisMaxFlips);
        
        % save some data
        
        KbQueueStop(k.kbDev);
        
        vbl = Screen('Flip', window);
        
        correct = response == stimTrialConds(3);
        
        data.data(t,:) = [b, t, stimTrialConds(1), stimTrialConds(2), stimTrialConds(3),...
            stimTrialConds(4), response, correct, rt, conf, confRT];
        data.prevResps{t} = prevResps;
        
        save(fullSaveFileName,'data','p');
        
    end
end

DrawFormattedText(window,'Thank you!\n\nThis part of the experiment is complete','center','center',p.black);
            
Screen('Flip',window); %display on screen

% this might be massive...
save(fullSaveFileName,'data','p','stimParamsS','stimParamsN','stimParamsC');



%--------------------------------------
%          3. Analyses
%---------------------------------------

% can add a summary of the data, if you like

% for example, 
% p(correct) in each condition
% mean(confidence) in each condition



WaitSecs(1);



%--------------------------------------
%          4. clean up
%---------------------------------------
       
clean_up
        

catch ME
    
    if exist('data','var')
        
        save(fullSaveFileName,'data','p','stimParamsS','stimParamsN','stimParamsC');
        
    end
    
    rethrow(ME)
    
end




end



