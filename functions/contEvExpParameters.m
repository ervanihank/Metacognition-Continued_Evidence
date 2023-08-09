function p = contEvExpParameters

p = struct;

p.dummyMode = 0;

% screen parameters

p.desiredRefreshRate = 100;

p.monitorSetUp = get(0, 'MonitorPositions');
p.screenNo = max(Screen('Screens'));
p.screenRefreshRate = Screen('FrameRate', p.screenNo);
p.flipInterval = 1/p.screenRefreshRate;

if p.dummyMode
    % if running on laptop
    p.flipInterval = 0.0167;
    p.desiredRefreshRate = 60;
end

% colour definitions
p.white = WhiteIndex(p.screenNo);
p.grey = p.white / 2;
p.black = BlackIndex(p.screenNo);

% screen res
% for the experiment computer
ScreenHypotenuse = 50.5;
pixHypotenuse = sqrt((p.monitorSetUp(3)^2)+(p.monitorSetUp(4)^2));
p.pixPerDeg = round(pixHypotenuse/ScreenHypotenuse);


% task parameters

% trials
% number of trials that are going to be repeated through the experiment
p.n_indiv_trials=90;
p.n_ranges = 3;
p.n_rep_range=p.n_indiv_trials/p.n_ranges;
% number of repetition of these individual trials
p.n_rep=1;
p.n_trials=p.n_indiv_trials*p.n_rep;

p.nBlocks = 6;

p.iti = 0.5;


p.directions = [0,180];
p.lifetime = 0.9;
p.win_size = 4; 

% maximum stim duration
p.maxFlips = 3*p.desiredRefreshRate;

p.p_noi = 0.05; 


% gabor parameters
% Gabor text Parameters
p.gabor_text_dim = round(p.pixPerDeg*0.75);
p.sigma = p.gabor_text_dim / 6;
p.contrast = 0.5 ;
p.aspectRatio = 1; 
% Spatial Frequency (Cycles Per Pixel)
p.num_cycles = 3;
p.freq = p.num_cycles / p.gabor_text_dim;
p.backgroundgabor_text = [0.5 0.5 0.5 0.5];
p.disable_norm = 1; 
p.pre_contrast_multiplier = 0.5; 


dim =(16/2)*p.pixPerDeg;
[x, y] = meshgrid(-dim:p.gabor_text_dim:dim, -dim:p.gabor_text_dim:dim);

%distance of each gabor_text from the center of the array
dist = sqrt(x.^2 + y.^2);

%Inner annulus
inner_dist = 1*p.pixPerDeg;
x(dist <= inner_dist) = nan;
y(dist <= inner_dist) = nan;

%Outer annulus
outer_dist = 8*p.pixPerDeg;
x(dist >= outer_dist) = nan;
y(dist >= outer_dist) = nan;

%Select only the finite values
x = x(isfinite(x));
y = y(isfinite(y));

% Center the annulus coordinates in the centre of the screen
if ~ p.dummyMode
    x_pos = x + 1920/2;
    y_pos = y + 1080/2;
else
    x_pos = x + 400;
    y_pos = y + 400;
end

% Count how many gabor_texts there are
p.n_gabors = numel(x_pos);

% Make the destination rectangles for all the gabor_texts in the array
base_rect = [0 0 p.gabor_text_dim p.gabor_text_dim];
all_rects = nan(4, p.n_gabors);
for i = 1:p.n_gabors
    all_rects(:, i) = CenterRectOnPointd(base_rect, x_pos(i), y_pos(i));
end

p.pos_each_gabor = all_rects;



% Global Drift speed
deg_per_sec = 3;
phasesPerCycle = 360;
p.drift_speed =  phasesPerCycle*deg_per_sec * (1/p.desiredRefreshRate);


% the ranges for each type of stimulus

p.rangeMaxs = [0.67,0.77,0.84];
p.rangeMin = 0.48;



