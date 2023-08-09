% function to perform one trial

% -----------------------------
%     INPUT
% -----------------------------

% window,           pointer to the window/screen object
% p,                parameters - loaded from contEvExpParameters
% k,                keyboard parameters - from main script
% gabor_tex,        procedural gabor
% orientations,     the orientations of each gabor
% phases            the phases of each gabor on each frame
% fixedDuration     whether the stimulus is presented for a fixed duration
% maxFlips          either the stimulus duration or the long maximum

% -----------------------------
%     OUTPUT
% -----------------------------

% response:    whether the participant pressed left or right
% rt:          when the button was pressed, relative to the start of the
%              stimulus

function [resp, rt, conf, confrt, prevResps] = presentTrialReplay(window, p, k, gabor_tex, orientations, phases, fixedDuration, maxFlips)


responded = 0;

KbReleaseWait(k.kbDev);
KbQueueFlush(k.kbDev);

prevResps = zeros(10,2); % up to 10 previous responses
pr_count = 1;

nextTrial = GetSecs + p.iti;

% in the iti - present the fixation point and do any computations

Screen('DrawDots', window, [p.center_x; p.center_y], 2, p.black, [], 2); %fixation point
Screen('DrawingFinished',window);
vbl = Screen('Flip', window);

properties=[NaN, p.freq, p.sigma, p.contrast, p.aspectRatio, 0, 0, 0];
propertiesMat = repmat(properties,p.n_gabors, 1, maxFlips);
propertiesMat(:,1,:)= phases;




while vbl<nextTrial
    
    Screen('DrawDots', window, [p.center_x; p.center_y], 2, p.black, [], 2); %fixation point
    Screen('DrawingFinished',window);
    vbl = Screen('Flip', window, vbl+(p.flipInterval*0.5)); % this tells the screen to flip asap
    
end

% set up the keyboard queue
KbQueueStart(k.kbDev);

fi = 1;

% now present the stimulus

startFlip = vbl+p.flipInterval;

doContinue = 1;

while doContinue
    
    
    % draw the stimulus first
    Screen('DrawTextures', window, gabor_tex, [], p.pos_each_gabor, orientations ,[], [], [], [], ...
        kPsychDontDoRotation, propertiesMat(:,:,fi)');
    Screen('DrawDots', window, [p.center_x; p.center_y], 2, p.black, [], 2); %fixation point
    Screen('DrawingFinished',window);
    vbl = Screen('Flip', window, vbl+(p.flipInterval*0.5));
    
    fi = fi+1;
    
    % check the keyboard
    
    [pressed, keys] = KbQueueCheck(k.kbDev);
            
    if pressed
        
        thisKey = min(find(keys,1,'first')); % gets the index of the first key in the list
        rt = min(keys(thisKey));
        
        switch thisKey
            case k.escapeKey
                fprintf('\n\nUser Quit')
                clean_up
                return
            case k.right % right arrow
                resp = 2; %leftwards is 1
                responded = 1;
            case k.left % left arrow
                resp = 1;
                responded = 1;
            otherwise % flush any other keys (not that they should be recorded)
                KbQueueFlush(k.kbDev);
        end
        
    end
    
    % check if we exit this routine
    
    if responded == 1
        
        if fixedDuration
            
            prevResps(pr_count,1) = resp;
            prevResps(pr_count,2) = rt-startFlip;
            pr_count = pr_count + 1;
            
            KbQueueFlush(k.kbDev);
            responded = 0;
        else
            doContinue = 0;
            
        end
            
    end
    
    if fi>maxFlips
        doContinue = 0;
    end
    
end

lastflip = vbl;


% if the participant still hasn't responded

while ~responded
    
    % present the response cue after a little break
    
    if GetSecs > (lastflip + p.iti/2)
        
        DrawFormattedText(window,'Enter final decision','center','center',p.black);
    else
        
        % keep presenting the fixation
        Screen('DrawDots', window, [p.center_x; p.center_y], 2, p.black, [], 2); %fixation point
        
    end
    
    Screen('DrawingFinished',window);
    vbl = Screen('Flip', window, vbl+(p.flipInterval*0.5));
    
    
    % check the keyboard
    
    [pressed, keys] = KbQueueCheck(k.kbDev);
            
    if pressed
        
        thisKey = min(find(keys,1,'first')); % gets the index of the first key in the list
        rt = min(keys(thisKey));
        
        switch thisKey
            case k.escapeKey
                fprintf('\n\nUser Quit')
                clean_up
                return
            case k.right % right arrow
                resp = 2; %leftwards is 1
                responded = 1;
            case k.left % left arrow
                resp = 1;
                responded = 1;
            otherwise % flush any other keys (not that they should be recorded)
                KbQueueFlush(k.kbDev);
        end
        
    end
    
end



% confidence decision

% wait a little bit
vbl = Screen('Flip', window, vbl+(p.flipInterval*0.5));

while GetSecs < (vbl + p.iti/2)
    % wait...
end

respondedConf = 0;

while ~respondedConf
    
        
    DrawFormattedText(window,'Confidence? 1 - 4','center','center',p.black);
    
    
    Screen('DrawingFinished',window);
    vbl = Screen('Flip', window, vbl+(p.flipInterval*0.5));
    
    
    % check the keyboard
    
    [pressed, keys] = KbQueueCheck(k.kbDev);
            
    if pressed
        
        thisKey = min(find(keys,1,'first')); % gets the index of the first key in the list
        confrt = min(keys(thisKey));
        
        switch thisKey
            case k.escapeKey
                fprintf('\n\nUser Quit')
                clean_up
                return
            case k.one 
                conf = 1; 
                respondedConf = 1;
            case k.two 
                conf = 2;
                respondedConf = 1;
            case k.three 
                conf = 3;
                respondedConf = 1;
            case k.four
                conf = 4;
                respondedConf = 1;
                
            otherwise % flush any other keys (not that they should be recorded)
                KbQueueFlush(k.kbDev);
        end
        
    end
    
end


% relative to the response
confrt = confrt - rt;

% and the reaction time relative to start
rt = rt - startFlip;

% trim the previous responses
prevResps = prevResps(1:pr_count,:);

end

    
    




    
    
    
    
    
    
    
    
    
    
    
    
