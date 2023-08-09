function practice(n_trials_prac)
% parameters about the fraction change
win_size=4; %signal fraction change every 4 frames
lifetime = 0.9; % the probability that any Gabor will change status
desiredRefreshRate=100;



% pixelperdeg=25;
% pixelperdeg=36;
% for the computer downstairs


% Open the screen
% [window, rect] = PsychImaging('OpenWindow', screenid, grey,[0,800;0,800]);
%full screen


%Screen('BlendFunction',window, 'GL_SRC_ALPHA','GL_ONE_MINUS_SRC_ALPHA')


% Refresh rate of the monitor






%% Number of trials
% n_trials=800;

%% Define gabor text Parameters

gabor_text_dim = round(pixelperdeg*0.75);

sigma = gabor_text_dim / 6;
contrast = 0.5 ;
aspectRatio = 1;
num_cycles = 3;
freq = num_cycles / gabor_text_dim;
% freq_one=num_cycles/gabor_text_dim_one;
backgroundgabor_text = [0.5 0.5 0.5 0.5];
disable_norm = 1;
pre_contrast_multiplier = 0.5;
gabor_text = CreateProceduralGabor(window, gabor_text_dim, gabor_text_dim,[],backgroundgabor_text , disable_norm, pre_contrast_multiplier);


% gabor_text_one=CreateProceduralGabor(window, gabor_text_dim_one, gabor_text_dim_one,[],backgroundgabor_text , disable_norm, pre_contrast_multiplier);



%% Position of gabor_texts
dim =(16/2)*pixelperdeg;
[x, y] = meshgrid(-dim:gabor_text_dim:dim, -dim:gabor_text_dim:dim); %create a matrix!!!!!

%distance of each gabor_text from the center of the array
dist = sqrt(x.^2 + y.^2);

%Inner annulus
inner_dist = 1*pixelperdeg;
x(dist <= inner_dist) = nan;
y(dist <= inner_dist) = nan;

%Outer annulus
outer_dist = 8*pixelperdeg;
x(dist >= outer_dist) = nan;
y(dist >= outer_dist) = nan;

%Select only the finite values
x = x(isfinite(x));
y = y(isfinite(y));

% Center the annulus coordinates in the centre of the screen
x_pos = x + center_x;
y_pos = y + center_y;

% Count how many gabor_texts there are
n_gabors = numel(x_pos);

% Make the destination rectangles for all the gabor_texts in the array
base_rect = [0 0 gabor_text_dim gabor_text_dim];
all_rects = nan(4, n_gabors);
for i = 1:n_gabors
    all_rects(:, i) = CenterRectOnPointd(base_rect, x_pos(i), y_pos(i));
end

%shuffle the positions
numgenerator=1:n_gabors;
pos_gabors=Shuffle(numgenerator);
pos_each_gabor=all_rects(1:4,pos_gabors);


%% Global Direction Parameters (global direction and speed)

%global direction (we have just 2 options,right/left, in the task)


% Global Drift speed
deg_per_sec = 3;
phasesPerCycle = 360;
drift_speed =  phasesPerCycle*deg_per_sec * (1/desiredRefreshRate);

%% Response Matrix
% 1 = response
% 2 = rt
% 3  = correctResponse

%% Setting Keyboard

kbDev = -1; % -1 for comp downstairs and for my own computer; 2 or 4 for the computer in lab
KbName('UnifyKeyNames');

%Keyboard Info
escapeKey = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');
space = KbName('space');

queueList = zeros(1,256);
queueList([escapeKey,leftKey,rightKey,space]) = 1;

KbQueueCreate(kbDev,queueList);

%% Timing, stimulus duration

% Interstimulus interval time in seconds and frames
isiTimeSecs = 0.5;
isiTimeFrames = round(isiTimeSecs / ifi);

% Numer of frames to wait before re-drawing
nextTrialStart = 0;
vbl = Screen('Flip', window);

%stimulus duration in flips
stimFlips = round(stimDur/ifi);

%%


%% parameters for matrix inversion
options = optimoptions('lsqlin', 'Display', 'off');
lb = [0, 0, 0, 0, 0, 0, 0, 0, 0];
ub = [1, 1, 1, 1, 1, 1, 1, 1, 1];

phases = zeros(n_gabors,1,stimFlips);


%% only one gabor 



pro_lst_practice1=(1:0.01:1);
pro_lst_practice2=(0.90:0.01:0.90);
pro_lst_practice3=(0.80:0.01:0.80);
pro_lst_practice4=(0.70:0.01:0.70);
pro_lst_practice5=(0.60:0.01:0.60);


n_range_prac=5; %number of range
rep_per_prolist_prac=2; %how many times one range will be presented
%create a list of range we want to test
list_prac=1:n_range_prac;
list_prac=repmat(list_prac,1,rep_per_prolist_prac);


%fraction of noise (stable for all trials)

stimFlips_prac=400;
directions_prac= [0,180,0,180,180,0,0,180,180,0]; %0 for left 180 for right
for trial=1:n_trials_prac
    
    a=list_prac(trial);
    if a==1
        pro_lst=pro_lst_practice1;
    elseif a==2
        pro_lst=pro_lst_practice2;
    elseif a==3
        pro_lst=pro_lst_practice3;
    elseif a==4
        pro_lst=pro_lst_practice4;
    elseif a==5
        pro_lst=pro_lst_practice5;
    end
    
    pro_nb=length(pro_lst);
       frac_transit = NaN(1, stimFlips_prac,n_trials_prac);
 gabor_state_mat = NaN(stimFlips_prac, n_gabors);
    %list of gabor driections
    direction_by_frame=NaN(stimFlips_prac, n_gabors);
    %list of phases
    phase_update= NaN(n_gabors,1,stimFlips_prac);
    
    orientations = randi([1,180],n_gabors,1); %orientation of each gabors will be stable within a trial
    direction=directions_prac(trial);
    contraDir = direction+180;
    
    propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],n_gabors, 1,stimFlips_prac);
    
    for tt = 1:stimFlips_prac
        if (mod(tt, win_size) == 1)
            
            pro_ind = randi(pro_nb);
            
            p_pro = pro_lst(pro_ind);
            p_noi=0;
            p_ant= 1-p_noi-p_pro;
            
            if (p_noi < 0)
                p_noi = 0;
                p_ant = 1 - p_pro;
            end
        end
        
        % -> build constraint matrix
        %    [a11, a12, a13, a21, a22, a23, a31, a32, a33]
        CC = [1, 1, 1, 0, 0, 0, 0, 0, 0; ...
            0, 0, 0, 1, 1, 1, 0, 0, 0; ...
            0, 0, 0, 0, 0, 0, 1, 1, 1; ...
            p_pro, 0, 0, p_ant, 0, 0, p_noi, 0, 0; ...
            0, p_pro, 0, 0, p_ant, 0, 0, p_noi, 0; ...
            0, 0, p_pro, 0, 0, p_ant, 0, 0, p_noi; ...
            p_pro, 0, 0, 0, p_ant, 0, 0, 0, p_noi];
        
        BB = [1, 1, 1, p_pro, p_ant, p_noi, lifetime]';
        
        AA = lsqlin(CC, BB, [], [], [], [], lb, ub, [], options);
        transit_mat = reshape(AA, 3, 3)';
        
        
        transit_count = 0;
        pro_count = 0;
        ant_count = 0;
        noi_count = 0;
        
        for ii = 1:n_gabors
            
            if (tt == 1)
                
                % -> initialize Gabor states and directions
                % -> choose state of Gabor: 1=signal, 2=anti-signal, 3=noise
                rand_val = rand;
                
                if (rand_val <= p_pro)
                    signal_val = 1;
                    dir=direction;
                elseif (rand_val <= (p_pro + p_ant))
                    signal_val = 2;
                    dir=contraDir;
                else
                    signal_val = 3;
                    dir=rand*360;
                end
                
                gabor_state_mat(1, ii) = signal_val;
                transit_detect = 1;
            else
                
                old_signal_val = gabor_state_mat(tt - 1, ii);
                rand_val = rand;
                
                switch old_signal_val
                    case 1
                        if (rand_val <= transit_mat(1,1))
                            signal_val = 1;
                            
                        elseif (rand_val <= (transit_mat(1,1) + transit_mat(1,2)))
                            signal_val = 2;
                            
                        else
                            signal_val = 3;
                            
                        end
                    case 2
                        if (rand_val <= transit_mat(2,1))
                            signal_val = 1;
                            
                        elseif (rand_val <= (transit_mat(2,1) + transit_mat(2,2)))
                            signal_val = 2;
                            
                        else
                            signal_val = 3;
                            
                            
                        end
                    case 3
                        if (rand_val <= transit_mat(3,1))
                            signal_val = 1;
                            
                        elseif (rand_val <= (transit_mat(3,1) + transit_mat(3,2)))
                            signal_val = 2;
                            
                        else
                            signal_val = 3;
                            
                        end
                    otherwise
                        fprintf('Error: %d state not defined\n', old_signal_val);
                end
                transit_detect = (signal_val ~= old_signal_val);
                transit_count = transit_count + (signal_val ~= old_signal_val);
                
            end
            if (transit_detect)
                switch signal_val
                    case 1
                        dir = direction;
                    case 2
                        dir = contraDir;
                    case 3
                        dir = rand*360;
                end
            else
                dir = direction_by_frame(tt-1, ii);
            end
            
            pro_count = pro_count + (signal_val == 1);
            ant_count = ant_count + (signal_val == 2);
            noi_count = noi_count + (signal_val == 3);
            
            gabor_state_mat(tt, ii) = signal_val;
            direction_by_frame(tt,ii)=dir;

        end
        
        phase_update(:,1,tt)=cosd(orientations'-direction_by_frame(tt,:)).* drift_speed;
        
        if (tt == 1)
            frac_transit(tt) = NaN;
        else
            frac_transit(tt) = transit_count / n_gabors;
        end

        if (tt==1)
            phases(:,1,tt) = randi([1,360],n_gabors,1,1); %starting phase
        else
            phases(:,1,tt)=phases(:,1,tt-1)+phase_update(:,1,tt);
        end
        
    end

    propertiesMat(:,1,:)= phases;
    
    
    KbReleaseWait(kbDev);
    KbQueueFlush(kbDev);
    KbQueueStart(kbDev);
    
    
    while GetSecs <= nextTrialStart
    end
    
    vbl = Screen('Flip', window);
    WaitSecs(0.4);
%     Screen('BlendFunction',window, 'GL_ONE','GL_ZERO');
    
    %
    tStart = GetSecs; %to count the RT
    fi = 1;
    pressed = 0;
    doContinue = 0;
    while ~doContinue
        while ~pressed
            [pressed, keys] = KbQueueCheck(kbDev);
            
            if fi < stimFlips_prac
                
                Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
                Screen('DrawTextures', window, gabor_text, [], pos_each_gabor, orientations' ,[], [], [], [], ...
                    kPsychDontDoRotation, propertiesMat(:,:,fi)');
                Screen('DrawDots', window, [center_x; center_y], 2, black, [], 2); %fixation point
                Screen('DrawingFinished',window);
                vbl = Screen('Flip', window);
                fi = fi+1;
 
            else
                blockMessage = ['What was the global direction?' '\n\n' ' \n\n ' 'Left? or Right?'];
                
                DrawFormattedText(window,blockMessage,'center','center',black);
                
                Screen('Flip',window); %display on screen

            end
        end
        
        thisKey = min(find(keys,1,'first'));
        keyTime = min(keys(thisKey));
        
        if thisKey == escapeKey
            fprintf('User Quit')
            ShowCursor;
            sca;
            return
        elseif thisKey == leftKey
            response = 1;
            if directions_prac(trial)==0
                correctResponse = 1;
                blockMessage = ['CORRECT ANSWER' '\n\n' 'It was moving leftward'];
                
                DrawFormattedText(window,blockMessage,'center','center',black);
                
                Screen('Flip',window); %display on screen
                WaitSecs(1);
                
            elseif directions_prac(trial)==180
                correctResponse=0;
                blockMessage = ['WRONG ANSWER' '\n\n' 'It was moving rightward'];
                
                DrawFormattedText(window,blockMessage,'center','center',black);
                
                Screen('Flip',window); %display on screen 
                 WaitSecs(1);
            end
            doContinue = 1;
        elseif thisKey == rightKey
            response = 2;
            if directions_prac(trial)==180
                correctResponse = 1;
                 blockMessage = ['CORRECT ANSWER' '\n\n' 'It was moving rightward'];
                
                DrawFormattedText(window,blockMessage,'center','center',black);
                
                Screen('Flip',window); %display on screen
                 WaitSecs(1);
                
            elseif directions_prac(trial)==0
                correctResponse=0;
                blockMessage = ['WRONG ANSWER' '\n\n' 'It was moving leftward'];
                
                DrawFormattedText(window,blockMessage,'center','center',black);
                
                Screen('Flip',window); %display on screen
                 WaitSecs(1);
                
            end
            doContinue = 1;
        else
            pressed = 0;
        end
        
    end
    
    KbQueueStop(kbDev);

    
    nextTrialStart = GetSecs+isiTimeSecs;
   
    
end
end