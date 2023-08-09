% function to generate stimulus parameters for the free task

% -----------------------------
%     INPUT
% -----------------------------


% n_gabors,             the number of gabors in the stimulus
% lifetime,             the lifetime for a direction
% win_size,             the window (flips) over which a stimulus probability 
%                       is sampled
% maxFlips,             maximum number of flips for the stimulus
% drift_speed,          the drift speed of the gabor (in the direction)
% p_noi,                the proportion of noise gabors
% directions,           the directions (usually [0,180], indexed by the
%                       conditionMat, position 3)
% rangeMin,             the minimum stimulus probability
% rangeMaxs,            that maximum for each range (indexed by the
%                       conditionMat, position 2)
% conditionMat          trial conditions, first column is the unique trial
%                       identifier, second column specifies the stimulus 
%                       probability range for that trial, third specifies 
%                       the direction

% -----------------------------
%     OUTPUT
% -----------------------------

% stimParams

% a structure containing (for each unique trial):

% trialID               unique identifier from condition matrix (typically 1:ntrials) 
% orientation           the orientation of each gabor
% phases                the phase of each gabor on each frame
% gabor_state_mat       whether the gabor is a signal (1), antisignal (2)
%                       or noise (3), for each frame
% direction_by_frame    the direction of each gabor on each frame
% dir                   the direction (index) the actual direction is
%                       directions(dir)
% val                   the random value for each gabor for each frame
% list_pro              the selected stimulus probabilities  for each frame
% rangeMax              the maximum stimulus probability on this trial
% thisevidence          the 'ideal evidence' of the stimulus on each frame


function [stimParams] = genStimParamsFree(n_gabors, lifetime, win_size, maxFlips, drift_speed,p_noi,directions, rangeMin, rangeMaxs,conditionMat)


ntrials = size(conditionMat,1);



stimParams = struct([]);
% ant_nb = length(ant_lst);
% -> parameters for matrix inversion
options = optimoptions('lsqlin', 'Display', 'off');
lb = [0, 0, 0, 0, 0, 0, 0, 0, 0];
ub = [1, 1, 1, 1, 1, 1, 1, 1, 1];


% go through the condition mat

for ti = 1:ntrials
    
    pro_lst = rangeMin:0.01:rangeMaxs(conditionMat(ti,2));
    pro_nb=length(pro_lst);
    
    direction = directions(conditionMat(ti,3));
    contraDir = direction+180;
    
    
    orientations = randi([1,180],1,n_gabors); %orientation of each gabors will be stable within a trial
    
    
    sigEv = cosd(orientations).*cosd(orientations);
    antiEv = cosd(orientations).*cosd(orientations)*-1;
    thisEvidence=zeros(maxFlips,1);
    
    
    % build the pro-list out here - to save time
    pro_inds = randi(pro_nb,ceil(maxFlips/win_size),1);
    pro_inds = repmat(pro_inds,1,win_size);
    pro_inds = reshape(pro_inds',maxFlips,1);
    
    list_pro = pro_lst(pro_inds);
    
    % pre-define signal vals
    val=rand(maxFlips,n_gabors);
    gabor_state_mat = ones(size(val))*3;
    gabor_state_mat(val <= list_pro(1)) = 1;
    gabor_state_mat(logical((val > list_pro(1)).*(val <= (1-list_pro(1)-p_noi)))) = 2;
    
    direction_by_frame = zeros(maxFlips,n_gabors);
    direction_by_frame(1,gabor_state_mat(1,:)==1) = direction;
    direction_by_frame(1,gabor_state_mat(1,:)==2) = contraDir;
    direction_by_frame(1,gabor_state_mat(1,:)==3) = rand(1,sum(gabor_state_mat(1,:)==3))*360;
    
    
    tet = 90-direction_by_frame(1,:);
    noiseEv = cosd(orientations).*cosd(tet'+orientations);
    
    thisEvidence(1,1) = (sum(sigEv(gabor_state_mat(1,:) == 1)) + sum(antiEv(gabor_state_mat(1,:) == 2)) + sum(noiseEv(gabor_state_mat(1,:) == 3)))/n_gabors;
    
    
    
    % now can start from the second flip
    
    
    for fi = 2:maxFlips
        
        p_pro=list_pro(fi);
        p_ant = 1-p_noi-p_pro;
        
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
        
        for gi = 1:n_gabors
            
            old_signal_val = gabor_state_mat(fi - 1, gi);
            rand_val=val(fi,gi);
            
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
                dir = direction_by_frame(fi-1, gi);
            end
            
            gabor_state_mat(fi, gi) = signal_val;
            direction_by_frame(fi,gi)=dir;
            
        end
        
        
        tet = 90-direction_by_frame(fi,:);
        noiseEv = cosd(orientations).*cosd(tet+orientations);
        
        thisEvidence(fi,1) = (sum(sigEv(gabor_state_mat(fi,:) == 1)) + sum(antiEv(gabor_state_mat(fi,:) == 2)) + sum(noiseEv(gabor_state_mat(fi,:) == 3)))/n_gabors;
        
    end
    
    phase_update = cosd(orientations-direction_by_frame).* drift_speed;
    
    startphases = randi([1,360],1,n_gabors);
    
    phases = startphases + cumsum(phase_update,1);
    
    
    stimParams(ti).trialID = conditionMat(ti,1);
    stimParams(ti).orientation = orientations;
    stimParams(ti).phases = round(phases,3);
    stimParams(ti).gabor_state_mat=gabor_state_mat;
    stimParams(ti).direction_by_frame=round(direction_by_frame,3);
    stimParams(ti).dir = conditionMat(ti,3);
    stimParams(ti).val = val;
    stimParams(ti).list_pro=list_pro;
    stimParams(ti).rangeMax = rangeMaxs(conditionMat(ti,2));
    stimParams(ti).thisevidence=round(thisEvidence,3);
    
end

end

            
            
            
    
    
    
    


