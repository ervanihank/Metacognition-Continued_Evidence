%% New Analysis Script
%Erva 
%25.03.22

close all
clear all

%% Load data

%Load free response block data first
dat_free= load([cd filesep 'data_contEvExp' filesep 'ENK' filesep 'ENK_contEvExp_freeTsk_S1_2022.mat']);
%Load replay block data
dat_replay= load([cd filesep 'data_contEvExp' filesep 'ENK_1' filesep 'ENK_1_contEvExp_replayTsk_S1_2022-03-24-23-21.mat']);
%stim params free
data_stim= load('stimParamsFree');
%% Required elemnets from Free Block
data_free=dat_free.data.data;
ifi=dat_free.p.flipInterval;
stimparams_free=dat_free.p;
duration_free=stimparams_free.maxFlips;
n_indiv_trials=stimparams_free.n_indiv_trials;
xvals_free = 1:duration_free;
n_trials_free=stimparams_free.n_trials;
drift_speed=stimparams_free.drift_speed;
n_rep_free=stimparams_free.n_rep;

%find which pre-defined stimulu represented in each trials in Free block
trial_id_free=data_free(:,3);

%Find maximum point of the range for each trials in the Free block
max_range_free=data_free(:,4);

%Response Free Block
resp_free=data_free(:,6);
correct_free=data_free(:,7);
rt_free=data_free(:,8);
rt_free(rt_free<0)=NaN;

%% Required elements from Replay Block

%% dectime

freeDec = zeros(n_indiv_trials,1);

for ti = 1:n_indiv_trials
    
    thesets = data_free(:,3) == ti;
    
    medrt = median(data_free(thesets,8));
    
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
freeDecFlips = round(freeDec/ifi);
lessFlips = round(lessDur/ifi);
moreFlips = round(moreDur/ifi);


decfram=freeDecFlips;

%%

stimparams_replay=dat_replay.p;
data_replay=dat_replay.data.data;

% load trial conditions in replay block
% 1=Free 2=Less 3=more- 4=More 5=More+
condition_replay=data_replay(:,6);
n_condition_type=length(unique(condition_replay));

%find which pre-defined stimulu represented in each trials in Replay block
trial_id_replay=data_replay(:,3);
n_trials_rep=stimparams_replay.n_trials;

n_gabors=stimparams_replay.n_gabors;

%Find maximum point of the range for each trials in the Replay block
max_range_replay=data_replay(:,4);

%responses Replay Block
prevResps_data=dat_replay.data.prevResps;
resp_replay=data_replay(:,7);
correct_replay=data_replay(:,8);
confi_replay=data_replay(:,10);
rt_replay=data_replay(:,9);
rt_replay(rt_replay<0)=NaN;

%initial responses in replay
for ii=1:n_trials_rep
resp_initial_replay(ii)=prevResps_data{ii,1}(1,1);
rt_initial(ii)=prevResps_data{ii,1}(1,2);
end

% number of range types in the replay block (it is same in the Free block)
range_type=unique(max_range_replay);
n_ranges=length(range_type);

%% Proportion Correct

% prop correct Free Block

for bb=1:n_ranges
    
    this_set=max_range_free==range_type(bb);
    p_correct_free(bb)=mean(correct_free(this_set));
    
    %calculate error for p correct
    n_trials_this_range=length(correct_free(this_set));
    sd_free(bb)=std(correct_free(this_set));
    err_free(bb)=sd_free(bb)/sqrt(n_trials_this_range);
    
end

% prop correct and confidence Replay Block

p_correct=nan(n_ranges,n_condition_type,1);
confidence=nan(n_ranges,n_condition_type,1);
sd_confi=[];
err_confi=[];
sd_pcorrect=[];
err_pcorrect=[];

for a=1:n_ranges
    
    for i=1:n_condition_type
        
        this_set=(condition_replay==i) &  (max_range_replay==range_type(a));
        n_trial_this_cond=length(confi_replay(this_set));
        p_correct(a,i)=mean(correct_replay(this_set));
        
        %error p_correct
        sd_pcorrect(a,i)=std(correct_replay(this_set));
        err_pcorrect(a,i)=sd_pcorrect(a,i)/sqrt(n_trial_this_cond);
        
        confidence(a,i)=mean(confi_replay(this_set));
        %error confidence
        sd_confi(a,i)=std(confi_replay(this_set));
        err_confi(a,i)=sd_confi(i)/sqrt(n_trial_this_cond);
        
    end
end

%% Evidence in the last frame
n_more_type=3;
n_rep_mores=2;
n_this_more=length(find(condition_replay==3));
last_evidence=nan(n_more_type,n_this_more,1);    

%last evidence in more
for b=1:90
    last_evidence_moreN(b)=dat_replay.stimParamsN(b).thisevidence(end);
end
   last_evidence(1,:)=repmat(last_evidence_moreN',n_rep_mores,1);
%last evidence in more+
for b=1:90
    last_evidence_moreS(b)=dat_replay.stimParamsS(b).thisevidence(end); 
end
last_evidence(2,:)=repmat(last_evidence_moreS',n_rep_mores,1);
     
%last evidence in more-
for b=1:90
    last_evidence_moreC(b)=dat_replay.stimParamsC(b).thisevidence(end);
end
last_evidence(3,:)=repmat(last_evidence_moreC',n_rep_mores,1);


stimparams=data_stim.stimParams;
last_evidence_free=nan(length(stimparams),1);

for ii=1:length(last_evidence_free)
       last_evidence_free(ii)=stimparams(ii).thisevidence(end);
end
        
last_evidence_free=repmat(last_evidence_free,n_rep_free,1);       

%% calculate how proportion correct change for each confidence rating

confi_type=4;
        
for dd=1:n_ranges
    for i=1:confi_type
        
        set_confi=(confi_replay==i) &  (max_range_replay==range_type(dd));
        n_this_confi(dd,i)=sum(set_confi);
        p_correct_confi(dd,i)=mean(correct_replay(set_confi));
        
    end
end

%for all trials independent of which range they are in
for i=1:confi_type
    
    set_confi=(confi_replay==i) ;
    n_this_confi_alltrials(i)=sum(set_confi);
    p_correct_confi_alltrials(i)=mean(correct_replay(set_confi));
    
end

%% rt in frees
n_trial_id=length(unique(trial_id_replay));

medRT_this_freeinfree=[];
medRT_this_freeinreplay=[];
medRT_this_initialinreplay=[];


for ff=1:n_trial_id
    
    set_this_freeinfree=trial_id_free==ff;
    medRT_this_freeinfree(ff)=nanmedian(rt_free(set_this_freeinfree));
    
    set_this_freeinreplay=(trial_id_replay==ff) & (condition_replay==1);
    medRT_this_freeinreplay(ff)=nanmedian(rt_replay(set_this_freeinreplay));
    
    %all first responses in the more conditions
    %if there is an initial resp take this 
    % if there is no initial resp take the normal resp
    cond_more= [condition_replay==3] + [condition_replay==4] + [condition_replay==5];
    set_this_initialinreplay=trial_id_replay==ff & cond_more==1;
    initials=rt_initial(set_this_initialinreplay);
    noinitial=find(initials==0);
    normalrts=rt_replay(set_this_initialinreplay);
    initials(noinitial)=normalrts(noinitial);
    medRT_this_initialinreplay(ff)=nanmedian(initials);

    
    %pure initial
    rt_pure_initial=rt_initial(set_this_initialinreplay);
    nanset=find(rt_pure_initial==0);
    rt_pure_initial(nanset)=nan;
    medRT_pure_initial(ff)=nanmedian(rt_pure_initial);
 
end  

%% change in response

resp_change_mores=[];

for i=3:5 %try to rewrite
    set_more=(condition_replay==i);
    n_mores(i)=sum(set_more);
    n_double_resp(i)=sum(rt_initial(set_more)~=0);
    n_single_resp(i)=sum(rt_initial(set_more)==0);
    
    %response change
    resp_initial_rep=resp_initial_replay(set_more);
    resp_rep=resp_replay(set_more);
    this_correct=correct_replay(set_more);
    fi=find(resp_initial_rep==0);
    resp_initial_rep(fi)=nan;
    resp_rep(fi)=nan;
    this_correct(fi)=nan;
    resp_initial_rep=resp_initial_rep(~isnan(resp_initial_rep));
    resp_rep=resp_rep(~isnan(resp_rep));
    this_correct=this_correct(~isnan(this_correct));

    resp_change_mores(i)=sum((resp_initial_rep)'~=resp_rep);
changing_set=(resp_initial_rep)'~=resp_rep;

p_correct_changing_after(i)=((sum(this_correct(changing_set)))/resp_change_mores(i))*100;
p_correct_changing_before(i)=100-p_correct_changing_after(i);
end
%OK

%% Plot

evidences=[0.10,0.15,0.20];

% define colors
mtlb_colors = [...
    0.0730, 0.5470, 0.6460;
    0.6400, 0.1730, 0.2950;
    0.9590, 0.6740, 0.6250;
    0.8940, 0.1850, 0.3560;
    0.0660, 0.3740, 0.1880;
    0.3010, 0.7550, 0.9330;
    0.2350, 0.0280, 0.1240;
    0.4350, 0.6780, 0.5840;
    0.6350, 0.0780, 0.1840;
    0.8000, 0,      0;
    0, 0, 0.5899;
    0.4567, 0.4498, 0.1212];

%% plot the evidence in the last frame 

figure
subplot(1,1,1)

hold on
%more
h1 = histogram(last_evidence(1,:),'Normalization','pdf');
h1.FaceAlpha = 0.3;
h1.NumBins=15;
%more+
h2 = histogram(last_evidence(2,:),'Normalization','pdf');
h2.FaceAlpha = 0.3;
h2.NumBins=15;
%more-
h3 = histogram(last_evidence(3,:),'Normalization','pdf');
h3.FaceAlpha = 0.3;
h3.NumBins=15;
%free
h4= histogram(last_evidence_free(:),'Normalization','pdf');
h4.FaceAlpha = 0.3;
h4.NumBins=15;

legend('More','More+','More-','Free')
xlabel('Evidence in the last frame','FontSize', 15)
ylabel('Density','FontSize', 15)

%% calculate gabor fraction status

%% Evidence and fraction changes
figure

% evidence change
subplot(1,3,1);
for i=1:n_condition_type

        find_ind=find(condition_replay==i & max_range_replay==range_type(1) & trial_id_replay==3);
        ind_evidence=find_ind(1);
        this_ID=trial_id_replay(ind_evidence);
        gabor_state_mat=stimparams(this_ID).gabor_state_mat;
        gabor_state_mat_size=size(gabor_state_mat);
        n_frames=gabor_state_mat_size(1);
        dec_time=decfram(this_ID);
        
        for tt=1:dec_time
            gabor_frac_status(1,tt) = (sum(gabor_state_mat(tt,:)==1)/n_gabors);
            gabor_frac_status(2,tt) = (sum(gabor_state_mat(tt,:)==2)/n_gabors);
            gabor_frac_status(3,tt) = (sum(gabor_state_mat(tt,:)==3)/n_gabors);
        end
        
        signal_fraction=gabor_frac_status(1,:);
        
        dec_time_in_rep=medRT_this_freeinreplay(this_ID)/ifi;
        dec_time_in_free=medRT_this_freeinfree(this_ID)/ifi;
        dec_time_in_initial=medRT_pure_initial(this_ID)/ifi;
        this_evidence=stimparams(this_ID).thisevidence;
        
        if i==1
            r=4;
            
        elseif i==2
            r=5;
        elseif i==3
            r=3;
        elseif i==4
            
            r=2;
        elseif i==5
            r=1;
        end
        
        if r==2
            
            evidence=this_evidence(1:80);
            
        elseif r==1
            less_dur=lessFlips(this_ID);
            evidence=this_evidence(1:less_dur);
            
        elseif r==3
            
            evidence=dat_replay.stimParamsC(this_ID).thisevidence;
            
        elseif r==5
            
            evidence=dat_replay.stimParamsN(this_ID).thisevidence;
            
        elseif r==4
            
            evidence=dat_replay.stimParamsS(this_ID).thisevidence;
            
        end
        
        
        this_duration= length(evidence);
        xvals = 1:this_duration;
        
        
        
        plot(xvals, evidence, '-' , 'Color', mtlb_colors(i,:),'LineWidth',2);
        hold on
        plot([dec_time,dec_time],[-0.5,0.5],'--','Color', mtlb_colors(10,:),'LineWidth',2);
        hold on
        
        plot([dec_time_in_rep,dec_time_in_rep],[-0.5,0.5],'--','Color', mtlb_colors(11,:),'LineWidth',2);
        hold on
        
        plot([dec_time_in_free,dec_time_in_free],[-0.5,0.5],'--','Color', mtlb_colors(7,:),'LineWidth',2);
        hold on
        plot([dec_time_in_initial,dec_time_in_initial],[-0.5,0.5],'--','Color', mtlb_colors(6,:),'LineWidth',2);
        hold on
end
title('evidence 0.10')
% xlabel('Frame')
ylabel('Evidence Replay')
xlim([0,50]);


% subplot(2,3,4)
% for i=1:n_condition_type
% 
%         find_ind=find(condition_replay==i & max_range_replay==range_type(1));
%         ind_evidence=find_ind(1);
%         this_ID=trial_id_replay(ind_evidence);
%         gabor_state_mat=stimparams(this_ID).gabor_state_mat;
%         for tt=1:n_frames
%             gabor_frac_status(1,tt) = (sum(gabor_state_mat(tt,:)==1)/n_gabors);
%             gabor_frac_status(2,tt) = (sum(gabor_state_mat(tt,:)==2)/n_gabors);
%             gabor_frac_status(3,tt) = (sum(gabor_state_mat(tt,:)==3)/n_gabors);
%         end
%         signal_fraction=gabor_frac_status(1,:);
%         
%         plot(xvals, signal_fraction, '-' , 'Color', mtlb_colors(i,:),'LineWidth',2);
%         hold on
%         
% end


% plot([dec_time,dec_time],[-1,1],'--','Color', mtlb_colors(10,:),'LineWidth',2);
% hold on
% plot([dec_time_in_rep,dec_time_in_rep],[-1,1],'--','Color', mtlb_colors(11,:),'LineWidth',2);
% hold on
% plot([dec_time_in_free,dec_time_in_free],[-1,1],'--','Color', mtlb_colors(7,:),'LineWidth',2);
% hold on
% plot([dec_time_in_initial,dec_time_in_initial],[-1,1],'--','Color', mtlb_colors(6,:),'LineWidth',2);
% hold on

% end

%     xlabel('Frame')
% ylabel('Fraction Replay')
% legend('more','more+','more-','freeinrep','less');
% ylim([0,1]);
% xlim([0,50]);

subplot(1,3,2);

for i=1:n_condition_type

        
        find_ind=find(condition_replay==i & max_range_replay==range_type(2) & trial_id_replay==31);
        ind_evidence=find_ind(1);
        this_ID=trial_id_replay(ind_evidence);
        gabor_state_mat=stimparams(this_ID).gabor_state_mat;
        gabor_state_mat_size=size(gabor_state_mat);
        n_frames=gabor_state_mat_size(1);
        dec_time=decfram(this_ID);
        
        for tt=1:dec_time
            gabor_frac_status(1,tt) = (sum(gabor_state_mat(tt,:)==1)/n_gabors);
            gabor_frac_status(2,tt) = (sum(gabor_state_mat(tt,:)==2)/n_gabors);
            gabor_frac_status(3,tt) = (sum(gabor_state_mat(tt,:)==3)/n_gabors);
        end
        
        signal_fraction=gabor_frac_status(1,:);
        
        dec_time_in_rep=medRT_this_freeinreplay(this_ID)/ifi;
        dec_time_in_free=medRT_this_freeinfree(this_ID)/ifi;
        dec_time_in_initial=medRT_pure_initial(this_ID)/ifi;
        this_evidence=stimparams(this_ID).thisevidence;
        
        if i==1
            r=4;
            
        elseif i==2
            r=5;
        elseif i==3
            r=3;
        elseif i==4
            
            r=2;
        elseif i==5
            r=1;
        end
        
        if r==2
            
            evidence=this_evidence(1:80);
            
        elseif r==1
            less_dur=lessFlips(this_ID);
            evidence=this_evidence(1:less_dur);
            
        elseif r==3
            
            evidence=dat_replay.stimParamsC(this_ID).thisevidence;
            
        elseif r==5
            
            evidence=dat_replay.stimParamsN(this_ID).thisevidence;
            
        elseif r==4
            
            evidence=dat_replay.stimParamsS(this_ID).thisevidence;
            
        end
        
        
        this_duration= length(evidence);
        xvals = 1:this_duration;
        
        
        
        plot(xvals, evidence, '-' , 'Color', mtlb_colors(i,:),'LineWidth',2);
        hold on
        plot([dec_time,dec_time],[-0.5,0.5],'--','Color', mtlb_colors(10,:),'LineWidth',2);
        hold on
        
        plot([dec_time_in_rep,dec_time_in_rep],[-0.5,0.5],'--','Color', mtlb_colors(11,:),'LineWidth',2);
        hold on
        
        plot([dec_time_in_free,dec_time_in_free],[-0.5,0.5],'--','Color', mtlb_colors(7,:),'LineWidth',2);
        hold on
        plot([dec_time_in_initial,dec_time_in_initial],[-0.5,0.5],'--','Color', mtlb_colors(6,:),'LineWidth',2);
        hold on
end
title('evidence 0.10')
% xlabel('Frame')
ylabel('Evidence Replay')
xlim([0,50]);

% subplot(2,3,5)
% 
% for i=1:n_condition_type
% 
%         find_ind=find(condition_replay==i & max_range_replay==range_type(2));
%         ind_evidence=find_ind(1);
%         this_ID=trial_id_replay(ind_evidence);
%         gabor_state_mat=stimparams(this_ID).gabor_state_mat;
%         for tt=1:n_frames
%             gabor_frac_status(1,tt) = (sum(gabor_state_mat(tt,:)==1)/n_gabors);
%             gabor_frac_status(2,tt) = (sum(gabor_state_mat(tt,:)==2)/n_gabors);
%             gabor_frac_status(3,tt) = (sum(gabor_state_mat(tt,:)==3)/n_gabors);
%         end
%         signal_fraction=gabor_frac_status(1,:);
%         
%         plot(xvals, signal_fraction, '-' , 'Color', mtlb_colors(i,:),'LineWidth',2);
%         hold on
%         
% end
% plot([dec_time,dec_time],[-1,1],'--','Color', mtlb_colors(10,:),'LineWidth',2);
% hold on
% plot([dec_time_in_rep,dec_time_in_rep],[-1,1],'--','Color', mtlb_colors(11,:),'LineWidth',2);
% hold on
% plot([dec_time_in_free,dec_time_in_free],[-1,1],'--','Color', mtlb_colors(7,:),'LineWidth',2);
% hold on
% plot([dec_time_in_initial,dec_time_in_initial],[-1,1],'--','Color', mtlb_colors(6,:),'LineWidth',2);
% hold on
% xlim([0,50]);
% 
% legend('freeinrep','less','more-','more','more+');
%   ylim([0,1]);
%   xlim([0,50]);


subplot(1,3,3);
for i=1:n_condition_type

        find_ind=find(condition_replay==i & max_range_replay==range_type(3) & trial_id_replay==61);
        ind_evidence=find_ind(1);
        this_ID=trial_id_replay(ind_evidence);
        gabor_state_mat=stimparams(this_ID).gabor_state_mat;
        gabor_state_mat_size=size(gabor_state_mat);
        n_frames=gabor_state_mat_size(1);
        dec_time=decfram(this_ID);
        
        for tt=1:dec_time
            gabor_frac_status(1,tt) = (sum(gabor_state_mat(tt,:)==1)/n_gabors);
            gabor_frac_status(2,tt) = (sum(gabor_state_mat(tt,:)==2)/n_gabors);
            gabor_frac_status(3,tt) = (sum(gabor_state_mat(tt,:)==3)/n_gabors);
        end
        
        signal_fraction=gabor_frac_status(1,:);
        
        dec_time_in_rep=medRT_this_freeinreplay(this_ID)/ifi;
        dec_time_in_free=medRT_this_freeinfree(this_ID)/ifi;
        dec_time_in_initial=medRT_pure_initial(this_ID)/ifi;
        this_evidence=stimparams(this_ID).thisevidence;
        
        if i==1
            r=4;
            
        elseif i==2
            r=5;
        elseif i==3
            r=3;
        elseif i==4
            
            r=2;
        elseif i==5
            r=1;
        end
        
        if r==2
            
            evidence=this_evidence(1:80);
            
        elseif r==1
            less_dur=lessFlips(this_ID);
            evidence=this_evidence(1:less_dur);
            
        elseif r==3
            
            evidence=dat_replay.stimParamsC(this_ID).thisevidence;
            
        elseif r==5
            
            evidence=dat_replay.stimParamsN(this_ID).thisevidence;
            
        elseif r==4
            
            evidence=dat_replay.stimParamsS(this_ID).thisevidence;
            
        end
        
        
        this_duration= length(evidence);
        xvals = 1:this_duration;
        
        
        
        plot(xvals, evidence, '-' , 'Color', mtlb_colors(i,:),'LineWidth',2);
        hold on
        plot([dec_time,dec_time],[-0.5,0.5],'--','Color', mtlb_colors(10,:),'LineWidth',2);
        hold on
        
        plot([dec_time_in_rep,dec_time_in_rep],[-0.5,0.5],'--','Color', mtlb_colors(11,:),'LineWidth',2);
        hold on
        
        plot([dec_time_in_free,dec_time_in_free],[-0.5,0.5],'--','Color', mtlb_colors(7,:),'LineWidth',2);
        hold on
        plot([dec_time_in_initial,dec_time_in_initial],[-0.5,0.5],'--','Color', mtlb_colors(6,:),'LineWidth',2);
        hold on
end
title('evidence 0.10')
% xlabel('Frame')
ylabel('Evidence Replay')
xlim([0,50]);

% subplot(2,3,6)
% for i=1:n_condition_type
% 
%         find_ind=find(condition_replay==i & max_range_replay==range_type(3));
%         ind_evidence=find_ind(1);
%         this_ID=trial_id_replay(ind_evidence);
%         gabor_state_mat=stimparams(this_ID).gabor_state_mat;
%         
%         for tt=1:n_frames
%             gabor_frac_status(1,tt) = (sum(gabor_state_mat(tt,:)==1)/n_gabors);
%             gabor_frac_status(2,tt) = (sum(gabor_state_mat(tt,:)==2)/n_gabors);
%             gabor_frac_status(3,tt) = (sum(gabor_state_mat(tt,:)==3)/n_gabors);
%         end
%         signal_fraction=gabor_frac_status(1,:);
%         
%         p1=plot(xvals, signal_fraction, '-' , 'Color', mtlb_colors(i,:),'LineWidth',2);
%         hold on
%         
% end
% 
% p2=plot([dec_time,dec_time],[-1,1],'--','Color', mtlb_colors(10,:),'LineWidth',2);
% hold on
% p3=plot([dec_time_in_rep,dec_time_in_rep],[-1,1],'--','Color', mtlb_colors(11,:),'LineWidth',2);
% hold on
% p4=plot([dec_time_in_free,dec_time_in_free],[-1,1],'--','Color', mtlb_colors(7,:),'LineWidth',2);
% hold on
% p5=plot([dec_time_in_initial,dec_time_in_initial],[-1,1],'--','Color', mtlb_colors(6,:),'LineWidth',2);
% hold on
% 
% xlim([0,50]);
% lh1=legend(p1,'freeinrep','less','more-','more','more+');
% lh2=legend([p2,p3,p4,p5],'dectime','median RTreplay','median RTfree','median RTinitialresp');
% ylim([0,1]);

%% Plot prop corrects and confidence
ranges=["evidence 0.10", "evidence 0.15", "evidence 0.20"];

figure
subplot

for i=1:n_ranges %for each range
    
    subplot(3,3,i)
    b=bar(p_correct(i,:));
      b.FaceColor = 'flat';
    hold on
    c=bar(6,p_correct_free(i));
     c.FaceColor = 'flat';
      C.CData(6,:) = mtlb_colors(6,:);

   for k = 1:5
    b.CData(k,:) = mtlb_colors(k,:);
   end
    
    hold on
    errorbar(1:5,p_correct(i,:),err_pcorrect(i,:),'Color', 'k','LineStyle','none')
 hold on
    errorbar(6,p_correct_free(i),err_free(i),'Color', 'k','LineStyle','none')
   
    hold on
    title(ranges(i))
    ylim([0.5, 1.1]);
 xlabel('condition')
  if i==1
    ylabel('Proportion Correct')
  end
  xticks([1 2 3 4 5 6])
   xticklabels({'FreeinReplay','Less','more-','more','more+','free'})
hold off
end

% plot averaged confidence for each conditions
hold on

for i=1:3
subplot(3,3,3+i)

b=bar(confidence(i,:));
b.FaceColor = 'flat';

   for k = 1:5
    b.CData(k,:) = mtlb_colors(k,:);
   end

hold on
errorbar(confidence(i,:),err_confi(i,:),'Color', 'k','LineStyle','none')

ylim([0, 4]);
xlabel('condition');
  if i==1
    ylabel('Average confidence')
  end
   xticklabels({'FreeinReplay','Less','more-','more','more+','free'})

hold off
end
%
    
for i=1:3
subplot(3,3,6+i)

b=bar(p_correct_confi(i,:));
b.FaceColor = 'flat';


for k = 1:4
    b.CData(k,:) = mtlb_colors(k+5,:);
end

hold on
xlabel('confidence level')
  if i==1
    ylabel('Proportion Correct')
  end
ylim([0.5, 1]);
xticklabels({'1','2','3','4'})
hold on
text(1:length(n_this_confi(i,:)),p_correct_confi(i,:),num2str(n_this_confi(i,:)'),'HorizontalAlignment','center','VerticalAlignment',"bottom"); 
hold off
end

%% metasensitivity
figure
subplot(1,1,1)

b=bar(p_correct_confi_alltrials);
b.FaceColor = 'flat';
for k = 1:4
    b.CData(k,:) = mtlb_colors(k+2,:);
end

hold on
xlabel('confidence level (all trials)')
ylabel('Proportion Correct (all trials)')

ylim([0.5, 1]);
xticklabels({'1','2','3','4'})
hold on
text(1:length(n_this_confi_alltrials),p_correct_confi_alltrials,num2str(n_this_confi_alltrials'),'HorizontalAlignment','center','VerticalAlignment',"bottom"); 
hold off

%% RT

median_rt_by_range=[];
mean_rt_by_range=[];
    
for ii=1:n_ranges
    
   thesets= max_range_free==range_type(ii);
%    thesets=thesets';
   
    %median rt for the given range
   median_rt_by_range(ii,1) = median(rt_free(logical(thesets.*(correct_free ==1))));
   median_rt_by_range(ii,2) = median(rt_free(logical(thesets.*(correct_free ==0))));
   median_rt_by_range(ii,3)= median(rt_free(thesets));
   
   %mean rt for the given range
   mean_rt_by_range(1,ii) = mean(rt_free(logical(thesets.*(correct_free ==1))));
   mean_rt_by_range(2,ii) = mean(rt_free(logical(thesets.*(correct_free ==0))));
   mean_rt_by_range(3,ii)= mean(rt_free(thesets));
    
   rt_free_by_range(ii,:)=rt_free(thesets);
   

end
    
%rt histogram for Free block for all trials
figure
subplot(2,2,1);

hold on
    

for i=1:n_ranges

h1=histogram(rt_free_by_range(i,:),'Normalization','pdf','FaceColor',mtlb_colors(i+2,:));
%h1.FaceAlpha = 0.1;
h1.NumBins=40;
hold on
end

xlabel('RT Distribution (for each range)')
ylabel('Density')
legend("evidence 0.10", "evidence 0.15", "evidence 0.20");
 
subplot(2,2,2)
hold on
h1=histogram(rt_free(:),'Normalization','pdf','FaceColor',mtlb_colors(6,:));
h1.NumBins=20;
xlabel('RT Distribution (all)')
ylabel('Density')

subplot(2,2,3);
hold on
bar( mean_rt_by_range(2,:));
bar( mean_rt_by_range(1,:));
legend('incorrect responses','correct responses')
xlabel('Evidence')
ylabel('mean rt ')
% ylim([0, 1.5]);
xticks([1 2 3])
xticklabels({'0.10', '0.15', '0.20'})

subplot(2,2,4)
hold on
bar( median_rt_by_range(:,2));
bar( median_rt_by_range(:,1));
legend('incorrect responses','correct responses')
xlabel('Evidence')
ylabel('median rt ')
% ylim([0, 1.5]);
xticks([1 2 3])
xticklabels({'0.10', '0.15', '0.20'})

figure
subplot(1,1,1)
hold on
bar(median_rt_by_range)
xlabel('evidence')
ylabel('median rt')
legend('Correct','Incorrect','All')
xticks([1 2 3])
xticklabels({'0.10', '0.15', '0.20'})
    
%% rt for each unique stimulus

figure
subplot(1,1,1)
hold on
all_trial_id=1:n_trial_id;
c1= mtlb_colors(1,:);
c2= mtlb_colors(2,:);
c3= mtlb_colors(3,:);
c4= mtlb_colors(4,:);
sz=25;
scatter(all_trial_id,medRT_this_freeinfree,sz,c1,'filled');
hold on
scatter(all_trial_id,medRT_this_freeinreplay,sz,c2,'filled');
hold on
scatter(all_trial_id,medRT_this_initialinreplay,sz,c3,'filled');
hold on
scatter(all_trial_id,medRT_pure_initial,sz,c4,'filled');

xlabel('unique stimulus number')
ylabel('median rt')
xlim([0, n_trial_id+1]);
legend('Free in Free','Free in Replay','Initial Resp in More','Pure Initial Resp More')

%% rt box plot for each range
data_this_freeinfree=[];
data_this_freeinreplay=[];
data_this_initialinreplay=[];
    
figure
hold on

for cc=1:n_ranges
    
      hold on
    set_this_freeinfree=max_range_free==range_type(cc);
    data_this_freeinfree=rt_free(set_this_freeinfree);  
    
    
    set_this_freeinreplay=(max_range_replay==range_type(cc)) & (condition_replay==1);
    data_this_freeinreplay=rt_replay(set_this_freeinreplay);
    
    cond_more= [condition_replay==3] + [condition_replay==4] + [condition_replay==5];
    set_this_initialinreplay=(max_range_replay==range_type(cc)) & cond_more==1;
    data_this_initialinreplay=rt_initial(set_this_initialinreplay);
    data_noinitial=find(data_this_initialinreplay==0);
    normalrtss=rt_replay(set_this_initialinreplay);
    data_this_initialinreplay(data_noinitial)=normalrtss(data_noinitial);
    
    
    %pure initial
    data_this_initialreplay_pure=rt_initial(set_this_initialinreplay);
    nanset=find(data_this_initialreplay_pure==0);
    data_this_initialreplay_pure(nanset)=nan;
    
     a=repmat({'freeinfree'},length(data_this_freeinfree),1);
     b=repmat({'freeinreplay'},length(data_this_freeinreplay),1);
     c=repmat({'initialinreplay'},length(data_this_initialinreplay),1);
     h=repmat({'pureinitial'},length(data_this_initialreplay_pure),1);
     x_axis=[a;b;c;h];
     ylabel('RT')
     data=[data_this_freeinfree;data_this_freeinreplay;data_this_initialinreplay';data_this_initialreplay_pure'];
     subplot(1,3,cc)

     boxplot(data,x_axis);
     title(ranges(cc))
     ylim([0 3]);
end

for zz=1:n_ranges
for rr=3:5

    more_conditions= condition_replay==rr & max_range_replay==range_type(zz);
    
    only_resp=rt_replay(more_conditions);
    
    this_initial=rt_initial(more_conditions);
    noansinitial=find(this_initial==0);
    normal_rep_rt=rt_replay(more_conditions);
    this_initial(noansinitial)=normal_rep_rt(noansinitial);
    this_initial=round(this_initial/ifi);
   
    
    
    this_stim_id=trial_id_replay(more_conditions);
    more_samps=moreFlips(this_stim_id);

    dec_frame=decfram(this_stim_id);
    
    for ii=1:length(dec_frame)
    if (dec_frame(ii)<=this_initial(ii)) && (this_initial(ii)<=more_samps(ii))
        dat_btw_dectimeANDmoretime(ii)=1;
    else
        dat_btw_dectimeANDmoretime(ii)=0;
    end
    end
    
for ii=1:length(dec_frame)
    if this_initial(ii)>more_samps(ii)
        dat_longer_than_moretime(ii)=1;
    else
        dat_longer_than_moretime(ii)=0;
    end
end
    
for ii=1:length(dec_frame)
    if this_initial(ii)<dec_frame(ii)
        before_dec(ii)=1;
    else
        before_dec(ii)=0;
    end
end

n_trials_btw_dectiemANDmoretime(zz,rr)=sum(dat_btw_dectimeANDmoretime);
n_trial_longer_than_moretime(zz,rr)=sum(dat_longer_than_moretime);
n_trial_before_dec(zz,rr)=sum(before_dec);

end
end


% plot the response changing in the more condition and the prop correct in
% that changes

figure
subplot(1,2,1)

bar([n_mores(3:5)',n_double_resp(3:5)',n_single_resp(3:5)',resp_change_mores(3:5)']);
xticks([1 2 3])
xticklabels({'more-','more','more+'})
xlabel('condition')
ylabel('count')
legend('number of more condition','number of double response','number of single response','number of response change')

subplot(1,2,2)
bar([p_correct_changing_before(3:5)',p_correct_changing_after(3:5)']);
hold on

hold on
xticks([1 2 3])
xticklabels({'more-','more','more+'})
xlabel('condition')
ylabel('prop correct')
legend('before response change','after response change')

%% 
figure
subplot(1,1,1)
hold on

data_freeinfree=rt_free;

cond_free=condition_replay==1;
data_freeinreplay=rt_replay(cond_free);

cond_more= [condition_replay==1] + [condition_replay==2] + [condition_replay==3];
initialinreplay_set=cond_more==1;
data_initialinreplay=rt_initial(initialinreplay_set);
data_noinitial=find(data_initialinreplay==0);
rts_last=rt_replay(initialinreplay_set);
data_initialinreplay(data_noinitial)=rts_last(data_noinitial);


data_initialreplay_pure=rt_initial(initialinreplay_set);
nanset=find(data_initialreplay_pure==0);
data_initialreplay_pure(nanset)=nan;

d=repmat({'freeinfree'},length(data_freeinfree),1);
e=repmat({'freeinreplay'},length(data_freeinreplay),1);
f=repmat({'initialinreplay'},length(data_initialinreplay),1);
g=repmat({'pureinitial'},length(data_initialreplay_pure),1);
x_axis=[d;e;f;g];
data_all_rts=[data_freeinfree;data_freeinreplay;data_initialinreplay';data_initialreplay_pure'];
ylabel('RT')
boxplot(data_all_rts,x_axis);

%%
figure

for i=1:3
subplot(1,3,i)


sz=25;
xax=[1 2 3 4];
yy=xax;

hold on
for aa=1:5

    
scatter(confidence(i,1),confidence(i,aa),sz,mtlb_colors(aa,:),'filled');
hold on

end
hold on
plot(xax,yy,'--', 'Color', mtlb_colors(1,:),'LineWidth',1);
xlim([1,4]);
ylim([1, 4]);
xlabel('confidence free in replay');
ylabel('Averaged confidence')
legend('FreeinReplay','Less','more-','more','more+')
title(num2str(evidences(i)))
end

%%
figure

for zz=1:n_ranges    

x=[];
for i=1:n_more_type
if i==1
    rr=3;
elseif i==2
    rr=4;
elseif i==4
    rr=5;
end
    a=zz+(3*(i-1));
    subplot(3,3,a);
    x(i,1)=n_trial_before_dec(zz,rr);
    
    x(i,2)=n_trials_btw_dectiemANDmoretime(zz,rr);
%     x(i,2)=n_trial_after_dec_before_passpoint(zz,i);
%     x(i,3)=n_trial_btw_passpointANDmoretime(zz,i);
    x(i,3)=n_trial_longer_than_moretime(zz,rr);
    
    pie(x(i,:));
    hold on
    
%     if a==1
%         ylabel('More-')
%     elseif a==4
%          ylabel('More')
%     elseif a==7
%          ylabel('More+')
%     end
    
    
    if a==1
        title(num2str(evidences(zz)))
    elseif a==2
        title(num2str(evidences(zz)))
    elseif a==3
        title(num2str(evidences(zz)))
    end

end
    
end
legend('Before dec time','Between dec time and more time','After more time');


%% just to show
figure
more_samps=moreFlips;
subplot(1,1,1);
for i=1:n_condition_type

        find_ind=find(condition_replay==i & max_range_replay==range_type(3) & trial_id_replay==61);
        ind_evidence=find_ind(1);
        this_ID=trial_id_replay(ind_evidence);
        gabor_state_mat=stimparams(this_ID).gabor_state_mat;
        gabor_state_mat_size=size(gabor_state_mat);
        n_frames=gabor_state_mat_size(1);
        dec_time=decfram(this_ID);
        
        for tt=1:dec_time
            gabor_frac_status(1,tt) = (sum(gabor_state_mat(tt,:)==1)/n_gabors);
            gabor_frac_status(2,tt) = (sum(gabor_state_mat(tt,:)==2)/n_gabors);
            gabor_frac_status(3,tt) = (sum(gabor_state_mat(tt,:)==3)/n_gabors);
        end
        
        signal_fraction=gabor_frac_status(1,:);
        
        dec_time_in_rep=medRT_this_freeinreplay(this_ID)/ifi;
        dec_time_in_free=medRT_this_freeinfree(this_ID)/ifi;
        dec_time_in_initial=medRT_pure_initial(this_ID)/ifi;
        this_evidence=stimparams(this_ID).thisevidence;
         more_time=more_samps(this_ID);
        
        if i==1
            r=4;
            
        elseif i==2
            r=5;
        elseif i==3
            r=3;
        elseif i==4
            
            r=2;
        elseif i==5
            r=1;
        end
        
        if r==2
            
            evidence=this_evidence(1:80);
            
        elseif r==1
            less_dur=lessFlips(this_ID);
            evidence=this_evidence(1:less_dur);
            
        elseif r==3
            
            evidence=dat_replay.stimParamsC(this_ID).thisevidence;
            
        elseif r==5
            
            evidence=dat_replay.stimParamsN(this_ID).thisevidence;
            
        elseif r==4
            
            evidence=dat_replay.stimParamsS(this_ID).thisevidence;
            
        end
        
        
        this_duration= length(evidence);
        xvals = 1:this_duration;


        plot(xvals, evidence, '-' , 'Color', mtlb_colors(i,:),'LineWidth',2);
        hold on

        plot([dec_time,dec_time],[-0.5,0.5],'--','Color', mtlb_colors(10,:),'LineWidth',3);
        hold on

        plot([dec_time_in_free,dec_time_in_free],[-0.5,0.5],'--','Color', mtlb_colors(7,:),'LineWidth',3);
        hold on
        plot([dec_time_in_initial,dec_time_in_initial],[-0.5,0.5],'--','Color', mtlb_colors(6,:),'LineWidth',3);
        hold on
        plot([more_time,more_time],[-0.5,0.5],'--','Color', mtlb_colors(12,:),'LineWidth',3);
        hold on


end

title('evidence 0.10')
% xlabel('Frame')
ylabel('Evidence Replay')
xlim([0,50]);
