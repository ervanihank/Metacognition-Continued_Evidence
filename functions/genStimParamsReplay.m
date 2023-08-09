% function to generate stimulus parameters for the replay task

% -----------------------------
%     INPUT
% -----------------------------



% -----------------------------
%     OUTPUT
% -----------------------------

% stimParams

% a structure containing (for each unique trial):
% orientation:      stimulus orientations
% phases:           the phase on each frame

function [stimParamsN, stimParamsS, stimParamsC] = genStimParamsReplay(stimParamsF, decFlips, moreFlips, n_gabors, rangeMin, p_noi, win_size,lifetime,drift_speed)

% More conditions
stimParamsN = struct([]);
stimParamsS = struct([]);
stimParamsC = struct([]);

% now saving less information in these matrices

directions = [0,180];

% ant_nb = length(ant_lst);
% -> parameters for matrix inversion
options = optimoptions('lsqlin', 'Display', 'off');
lb = [0, 0, 0, 0, 0, 0, 0, 0, 0];
ub = [1, 1, 1, 1, 1, 1, 1, 1, 1];


nts = length(stimParamsF);

for ti = 1:nts
    
    
    % redefine
    maxFlips = decFlips(ti) + moreFlips(ti);
    
    phases = zeros(maxFlips,n_gabors);
    gabor_state_mat = zeros(maxFlips,n_gabors);
    direction_by_frame = zeros(maxFlips,n_gabors);
    val=rand(maxFlips,n_gabors);
    
    
    
    
    
    thisEvidence = zeros(maxFlips,1);
    
    
    
    % previous trial up to decision
    orientations = stimParamsF(ti).orientation;
    direction = directions(stimParamsF(ti).dir);
    contraDir = direction + 180;
    
    phases(1:decFlips(ti),:) = stimParamsF(ti).phases(1:decFlips(ti),:);
    direction_by_frame(1:decFlips(ti),:) = stimParamsF(ti).direction_by_frame(1:decFlips(ti),:);
    gabor_state_mat(1:decFlips(ti),:) = stimParamsF(ti).gabor_state_mat(1:decFlips(ti),:);
    
    val(1:decFlips(ti),:) = stimParamsF(ti).val(1:decFlips(ti),:);
    thisEvidence(1:decFlips(ti),:) = stimParamsF(ti).thisevidence(1:decFlips(ti),:);
    
    sigEv = cosd(orientations).*cosd(orientations);
    antiEv = cosd(orientations).*cosd(orientations)*-1;
    
    
    % adding more - neutral condition
    
    pro_lst = rangeMin:0.01:stimParamsF(ti).rangeMax;
    pro_nb = length(pro_lst);
    
    pro_inds = randi(pro_nb,ceil(maxFlips/win_size),1);
    pro_inds = repmat(pro_inds,1,win_size);
    pro_inds = reshape(pro_inds',ceil(maxFlips/win_size)*win_size,1);
    
    list_pro = pro_lst(pro_inds);
    list_pro(1:decFlips(ti)) = stimParamsF(ti).list_pro(1:decFlips(ti));
    
    % just go through the additional frames
    
    for fi = (decFlips(ti)+1):maxFlips
        
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
        
        thisEvidence(fi,1) = sum(sigEv(gabor_state_mat(fi,:) == 1)) + sum(antiEv(gabor_state_mat(fi,:) == 2)) + sum(noiseEv(gabor_state_mat(fi,:) == 3))/n_gabors;
        
    end
    
    phase_update = cosd(orientations-direction_by_frame).* drift_speed;
    
    startphases = phases(1,:) - phase_update(1,:);
    
    phases = startphases + cumsum(phase_update,1);
    
    stimParamsN(ti).trialID = stimParamsF(ti).trialID;
    stimParamsN(ti).orientation = orientations;
    stimParamsN(ti).phases = round(phases,3);
    %stimParamsN(ti).gabor_state_mat=gabor_state_mat;
    stimParamsN(ti).direction_by_frame=round(direction_by_frame,3);
    %stimParamsN(ti).dir = stimParamsF(ti).dir;
    %stimParamsN(ti).val = val;
    %stimParamsN(ti).list_pro=list_pro;
    %stimParamsN(ti).rangeMax = stimParamsF(ti).rangeMax;
    stimParamsN(ti).thisevidence=round(thisEvidence,3);
    
    
    
    
    
    % adding more - supporting evidence condition
    
    % new random numbers
    val=rand(maxFlips,n_gabors);
    val(1:decFlips(ti),:) = stimParamsF(ti).val(1:decFlips(ti),:);
    
    
    pro_lst = stimParamsF(ti).rangeMax:0.01:0.90;
    pro_nb = length(pro_lst);
    
    pro_inds = randi(pro_nb,ceil(maxFlips/win_size),1);
    pro_inds = repmat(pro_inds,1,win_size);
    pro_inds = reshape(pro_inds',ceil(maxFlips/win_size)*win_size,1);
    
    list_pro = pro_lst(pro_inds);
    list_pro(1:decFlips(ti)) = stimParamsF(ti).list_pro(1:decFlips(ti));
    
    % just go through the additional frames
    
    for fi = (decFlips(ti)+1):maxFlips
        
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
    
    startphases = phases(1,:) - phase_update(1,:);
    
    phases = startphases + cumsum(phase_update,1);
    
    stimParamsS(ti).trialID = stimParamsF(ti).trialID;
    stimParamsS(ti).orientation = orientations;
    stimParamsS(ti).phases = round(phases,3);
    %stimParamsS(ti).gabor_state_mat=gabor_state_mat;
    stimParamsS(ti).direction_by_frame=round(direction_by_frame,3);
    %stimParamsS(ti).dir = stimParamsF(ti).dir;
    %stimParamsS(ti).val = val;
    %stimParamsS(ti).list_pro=list_pro;
    %stimParamsS(ti).rangeMax = stimParamsF(ti).rangeMax;
    stimParamsS(ti).thisevidence=round(thisEvidence,3);
    
    
    
    
    
    
    
    
    
    
    
    % adding more - contrary evidence condition
    
    % new random numbers
    val=rand(maxFlips,n_gabors);
    val(1:decFlips(ti),:) = stimParamsF(ti).val(1:decFlips(ti),:);
    
    
    
    pro_lst = 0.10:0.01:rangeMin;
    pro_nb = length(pro_lst);
    
    pro_inds = randi(pro_nb,ceil(maxFlips/win_size),1);
    pro_inds = repmat(pro_inds,1,win_size);
    pro_inds = reshape(pro_inds',ceil(maxFlips/win_size)*win_size,1);
    
    list_pro = pro_lst(pro_inds);
    list_pro(1:decFlips(ti)) = stimParamsF(ti).list_pro(1:decFlips(ti));
    
    % just go through the additional frames
    
    for fi = (decFlips(ti)+1):maxFlips
        
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
    
    startphases = phases(1,:) - phase_update(1,:);
    
    phases = startphases + cumsum(phase_update,1);
    
    stimParamsC(ti).trialID = stimParamsF(ti).trialID;
    stimParamsC(ti).orientation = orientations;
    stimParamsC(ti).phases = round(phases,3);
    %stimParamsC(ti).gabor_state_mat=gabor_state_mat;
    stimParamsC(ti).direction_by_frame=round(direction_by_frame,3);
    %stimParamsC(ti).dir = stimParamsF(ti).dir;
    %stimParamsC(ti).val = val;
    %stimParamsC(ti).list_pro=list_pro;
    %stimParamsC(ti).rangeMax = stimParamsF(ti).rangeMax;
    stimParamsC(ti).thisevidence=round(thisEvidence,3);
    
    
end
    
    
    

