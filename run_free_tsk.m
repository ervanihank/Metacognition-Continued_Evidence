% script to run Free Task with drifting gabors stimuli
Screen('Preference', 'SkipSyncTests', 1);
% each trial the observer is shown an arrage of gabors, whose phase is
% updated to appear as if drifting in a particular direction
% the observer must decide as quickly as possible if the gabors are
% drifting leftward or rightward. The stimuli are displayed until the
% response is entered.

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

% none
% stimParams are loaded from the saved file
%

%-----------------------------------
%       Output
%-----------------------------------

% p:        parameters
% data:     data 



%-----------------------------------
%       1. Set up
%------------------------------------

try

commandwindow;

experimentName = 'contEvExp';
taskName = 'freeTsk';

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


% Check if the data folder exists already
data_folder = ['data_' experimentName];

[ t, newDir ] = isDirectory(data_folder);

if t
    disp(['Main data save directory set to: ', data_folder, '']);
else
    disp(['\nOops! Cannot find main data save directory: ', data_folder, '\n']);
    make_folder = input('Would you like to make this directory? y or n: ', 's' );
    if make_folder == 'y'
        mkdir(data_folder);
    else
        fprintf('\nExiting script now...\n\n')
        return;
    end
end


% Participant Number

% Display subjects already used
fprintf('\nThese are the existing subject data directories:\n\n');
ls(data_folder);

% Get the subject number
subject = upper( input('Please enter the participant''s ID code [Test]:     ', 's') );

if isempty(subject)
    subject = 'Test';
    disp(['No subject entered, defaulting to ', subject, ' directory.']);
end

% Data storage

dataSavePath = [data_folder filesep subject];

% check for participant folder
[ t, newDir ] = isDirectory(dataSavePath);
if ~t
    createNewDir = upper( input(['\nDo you want to create a new data save directory for ', subject, ' (y/n)? '], 's'));
    if createNewDir == 'Y'
        mkdir(dataSavePath);
        fprintf('\nOkay, new directory created:\n\n');
        ls(data_folder);
    else
        fprintf('\nOkay, exiting the script now...\n\n');
        return;
    end
end

% previous data files
dataFiles = ls(dataSavePath);
if ~isempty(dataFiles)
    fprintf(['\n\nHere are the existing data files for ', subject, ':\n\n']);
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

if validFiles == 0
    fprintf('\nThere are currently no valid run data filenames for this subject...\n');
    session = 1;
    
else
    fprintf('\n');
    for f = 1:validFiles
        fprintf(['(', int2str(f), ')  ', runFiles{f}, '\n']);
    end
    session = ( input('Please enter the session number:     ') );
end


sessiontype = ['S' num2str(session)];

fprintf('\n-------------------------------\n');
fprintf([' THIS RUN IS SESSION ' int2str(session) '\n']);
fprintf('-------------------------------\n\n');


% Get date and time for trial
time_stamp = datestr(now, 31);       % get current date and time as string
time_stamp = time_stamp(1:end-3);    % chop seconds off
time_stamp(end-[2 5]) = '-';         % replace space and colon with dashes

% Name files
dataSaveFileName = strcat(subject, '_', experimentName, '_', taskName, '_', sessiontype , '_', time_stamp);
fullSaveFileName = [dataSavePath filesep dataSaveFileName '.mat'];  % final filename with full path
disp(['New data filename is: ', fullSaveFileName]);    
    


% parameters and conditions 

% load parameters from function script

p = contEvExpParameters;

p.experimentName = experimentName;
p.participant = subject;
p.session = session;

% set up random seed to be able to regenerate experiment
p.startTime = time_stamp;
p.dataFile = fullSaveFileName;
p.taskName = taskName;

p.initRandSeed = GetSecs;
% rng('default');
rng(p.initRandSeed);


% conditions

% this is linked to the set up in genStimParamsFree

% the conditions are:
% which stimulus
% which range is that stimulus in
% which direction is it going

stimIdList = 1:p.n_indiv_trials;
rangeList = sort(repmat(1:p.n_ranges,1,p.n_indiv_trials/p.n_ranges));
dirList = repmat(sort(repmat(1:2,1,p.n_indiv_trials/p.n_ranges/2)),1,p.n_ranges);

% now each stimulus is linked to a range and a direction (and there are
% equal directions for each range - to avoid biases)

p.conditionMat = [stimIdList',rangeList',dirList'];

% the participant is shown this set of trials p.n_rep times
p.trialConditions = repmat(p.conditionMat,p.n_rep,1);

% shuffle the order
shuffInds = randperm(p.n_trials);
p.trialConditions = p.trialConditions(shuffInds,:);

% % generating evidence for n_repeated_trials to be repeated for n_rep times    
% [stimParams] = genStimParamsFree(p.n_gabors, p.lifetime, p.win_size,...
%     p.maxFlips, p.drift_speed, p.p_noi, p.directions, p.rangeMin,...
%     p.rangeMaxs, p.conditionMat);

% load stim params from file
% (to re-make the stimulus parameters, use the function above)

load([cd filesep 'stimParamsFree.mat'])


% data storage

data.headers = {'block'; 'Trial'; 'stimID'; 'range'; 'trueDir'; 'Response'; 'Correct'; 'RT'};
data.data = zeros(p.n_trials, length(data.headers));



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
k.escapeKey = KbName('ESCAPE');
k.right = KbName('rightArrow');
k.left = KbName('leftArrow');
k.space = KbName('space');

ListenChar(-1);

% List for KbQueue
k.queueList = zeros(1,256);
k.queueList([k.escapeKey, k.right, k.left, k.space]) = 1;

KbQueueCreate(k.kbDev,k.queueList);


% stimuli

% create the procedural gabor
gabor_tex = CreateProceduralGabor(window, p.gabor_text_dim, p.gabor_text_dim,[],p.backgroundgabor_text , p.disable_norm, p.pre_contrast_multiplier);


%DrawFormattedText(window,'Loading experiment...','center','center',p.black);
    
%Screen('Flip',window); %display on screen





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
        
        theseOrientations = stimParams(stimTrialConds(1)).orientation';
        thesePhases = stimParams(stimTrialConds(1)).phases';
        
        
        % present the stimulus and get the response
        
        [response, rt] = presentTrialFree(window, p, k, gabor_tex, theseOrientations, thesePhases);
        
        % generate fake data
        %response = randi(2,1);
        %rt = abs(randn(1)+0.2);
        
        
        % save some data
        
        KbQueueStop(k.kbDev);
        
        vbl = Screen('Flip', window);
        
        correct = response == stimTrialConds(3);
        
        data.data(t,:) = [b, t, stimTrialConds(1), stimTrialConds(2), stimTrialConds(3), response, correct, rt];
        save(fullSaveFileName,'data','p');
        
    end
end

DrawFormattedText(window,'Thank you!\n\nThis part of the experiment is complete','center','center',p.black);
            
Screen('Flip',window); %display on screen

save(fullSaveFileName,'data','p');


%--------------------------------------
%          3. Analyses
%---------------------------------------

% can add a summary of the data

% for example, 
% p(correct) in each range
% rt in each range

pcorr_range = zeros(1,p.n_ranges);
nbins = 20;
rt_dist = zeros(p.n_ranges,nbins);
rt_x = zeros(p.n_ranges,nbins);

for ri = 1:p.n_ranges
    thesets = data.data(:,3)==ri;
    
    pcorr_range(ri) = mean(data.data(thesets,7));
    
    [di,edg] = histcounts(data.data(thesets,8),nbins);
    rt_dist(ri,:) = di;
    rt_x(ri,:) = edg(1:nbins)+((edg(2)-edg(1))/2);
end


figure;
subplot(2,1,1)
bar(pcorr_range)
axis([0.5 p.n_ranges+0.5 0.5 1])
xlabel('Range condition')
ylabel('Proportion correct')

subplot(2,1,2)
hold on
for ri = 1:p.n_ranges
    plot(rt_x(ri,:),rt_dist(ri,:),'LineWidth',2)
end
xlabel('Reaction time')
ylabel('Count')



WaitSecs(1);



%--------------------------------------
%          4. clean up
%---------------------------------------
       
clean_up
        

catch ME
    
    if exist('data','var')
        
        save(fullSaveFileName,'data','p');
        
    end
    
    rethrow(ME)
    
end
        
        
        
        



















