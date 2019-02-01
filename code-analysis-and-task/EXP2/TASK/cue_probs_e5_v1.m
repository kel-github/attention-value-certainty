%%%%%%%%%%%%% K. Garner - Sept 2016, University of Birmingham
%%%%%%%%%%%%% task is a relative value cueing task where participants respond to a cued
%%%%%%%%%%%%% target (.8/.2, .6/.4), irrespective of
%%%%%%%%%%%%% the reward that target's placeholder represents.
%%%%%%%%%%%%% Cue-target SOA = 50 ms (Lanthier et al, 2015).
%%%%%%%%%%%%% differs from exp 1 in that there is now one condition where
%%%%%%%%%%%%% an exponential decay is applied to the reward value upon
%%%%%%%%%%%%% target onset (decay vs static reward conditions)


% _______________________________________________________________________________________________________________
% DEPENDENCIES
% coded on Matlab 2013a using PTB v3.0.12 
% run on Stone SOFREP-144
% with 23-inch Asus VG278HE monitor (1920 x 1080 pixels, 60-Hz refresh)

%%%%% clear old stuff
clear all
clear mex

% initialise mex files etc
KbCheck;
KbName('UnifyKeyNames');
GetSecs;
AssertOpenGL
Screen('Preference', 'SkipSyncTests', 0);
%PsychDebugWindowConfiguration;
dbug = input('Debug Mode? (1 or 0) ');

if any(dbug)
    
    sess.sub_num = 999;
    sess.blocks = [1 4 2 3];
    sess.response_order = 1;
    sess.shapes = [1 2];
else
    
    sess.sub_num = input('Subject Number? ');
    load(sprintf('sub_info_%d', sess.sub_num));
    sess.blocks = sub_info(4:7); % 1 = .4/.6 static, 2 = .4/.6 decay, 3 = .2/.8 static, 4 = .2/.8 decay
    sess.response_order = sub_info(3);
    if sub_info(2) == 1
        sess.shapes = [1 2];
    else
        sess.shapes = [2 1];
    end
end
sess.eye_on = input('Eyetracking? ');

%%%%%%%%% set time cheats
expand_time = 1;

% set randomisation seed based on sub/sess number
sess.session = 1;
r_num = [num2str(sess.sub_num) num2str(sess.session)];
r_num = str2double(r_num);
rand('state',r_num);
randstate = rand('state');

% assign save mat
save_mat = ['KG_exp3_cueprob_v5_1_' num2str(sess.sub_num) '_s' num2str(sess.session)];

parent = cd;
% local = [num2str(sess.sub_num) 's' num2str(sess.session)];
% writeDir(parent,local);
% full_f_name = fullfile(parent,local,save_mat);

%%%%%%%%% set up screen and draw stimuli
screen_nums = Screen('Screens');
screen_nums = max(screen_nums);
% get black and white
white = WhiteIndex(screen_nums);
black = BlackIndex(screen_nums);
screen_grey = 118;
[w, rect] = Screen('OpenWindow', screen_nums ); % ,...
[ screen_x_pix, screen_y_pix ] = Screen( 'WindowSize', w );
Screen(w, 'Flip');
Screen(w, 'FillRect', screen_grey, rect);
[x_center, y_center] = RectCenter(rect);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
ifi = Screen('GetFlipInterval', w);
% Retreive the maximum priority number
topPriorityLevel = MaxPriority(w);

%%%%%%% draw stimuli
%%%%%% square vs circle (blue) - define rects and colours
tgt_size = 50; % start with 30 and amend accordingly
center_offset = 150;
tgt_base_rects = [0 0 tgt_size tgt_size];
tgt_x_coords = [x_center - center_offset, x_center + center_offset];
num_tgt_locs = 2;
tgt_rect_coords = zeros(4, num_tgt_locs);
for i = 1:num_tgt_locs
    tgt_rect_coords(:, i) = CenterRectOnPointd(tgt_base_rects, tgt_x_coords(i), y_center);
end
tmp = [86.7 74.562 80.631; 119.85 86.2321 1.1985]; % check colour luminance
task.shape_col = [tmp(sess.shapes(1),:); tmp(sess.shapes(2),:)]; % now call 1 for no reward, 2 for high
oval_width = 8;
neut_col = [155 155 155];
stooge_col = [ 120 120 120 ; 120 120 120 ];

%%%%%%%% now draw the targets (letters H or N)
target_col = [ 90 90 90];
target_center_left = x_center - center_offset; %%%%%% place to center target on
target_center_right = x_center + center_offset;
tgts = {'H','N'};
dsts = {'K','Z'};
x_offset = 8;
y_offset = 15;
tgt_y = y_center-y_offset;
dst_y = y_center-y_offset;


%%%%% draw two fixations- consisting of two straight lines or two arrows, one pointing left, one
%%%%% right
alert_fix_pix = 6;
alert_fix_off = 2;
fix_a_x = [-alert_fix_off -(alert_fix_off+alert_fix_pix) -(alert_fix_off+alert_fix_pix) -alert_fix_off alert_fix_off alert_fix_off+alert_fix_pix alert_fix_off+alert_fix_pix alert_fix_off];
fix_a_y = [-alert_fix_pix 0 0 alert_fix_pix -alert_fix_pix 0 0 alert_fix_pix];
fix_a = [fix_a_x; fix_a_y];
fix_b_x = [-((alert_fix_off+alert_fix_pix)/2) -((alert_fix_off+alert_fix_pix)/2) (alert_fix_off+alert_fix_pix)/2 (alert_fix_off+alert_fix_pix)/2];
fix_b_y = [-alert_fix_pix alert_fix_pix -alert_fix_pix alert_fix_pix];
fix_b = [fix_b_x; fix_b_y];
fix_all = {fix_a; fix_b};
fix_allocate = [2 2]; % this refers to earlier versions where fixation changed
fix_all = {fix_all{fix_allocate(1)} fix_all{fix_allocate(2)}}; %%%%%%%%%%% randomising first and second fix across participants

base_cols_a = [155 155 155; 155 155 155; 155 155 155; 155 155 155; 155 155 155;155 155 155; 155 155 155; 155 155 155]';
base_cols_b = [155 155 155; 155 155 155; 155 155 155; 155 155 155]';
alert_base_cols_all = {base_cols_a base_cols_b};
alert_base_cols_all = {alert_base_cols_all{fix_allocate(1)} alert_base_cols_all{fix_allocate(2)}};
alert_fix_cols_a = {[100 100 100;100 100 100; 100 100 100; 100 100 100; 155 155 155; 155 155 155; 155 155 155; 155 155 155]',[155 155 155; 155 155 155; 155 155 155; 155 155 155; 100 100 100; 100 100 100; 100 100 100; 100 100 100]'};
alert_fix_cols_b = {[100 100 100; 100 100 100; 155 155 155; 155 155 155]',[155 155 155; 155 155 155; 100 100 100; 100 100 100]'};
alert_fix_cols_all = {alert_fix_cols_a alert_fix_cols_b};
alert_fix_cols_all = {alert_fix_cols_all{fix_allocate(1)} alert_fix_cols_all{fix_allocate(2)}};
alert_fix_width = 4;
reward_fix_width = 8;
alert_fix_check_coords = CenterRectOnPointd([0 0 100 100], x_center, y_center); %%%% for eyetracking

%%%%%% drawing stuff for measure of perceived value distance screen
vas_coords = [-500 500; 0 0];
vas_width_pix = 6;

vas_tgt_y_coords = [y_center - 150, y_center - 150];
num_tgt_locs = 2;
vas_tgt_rect_coords = zeros(4, num_tgt_locs);
vas_tgt_rect_coords(:, 1) = CenterRectOnPointd(tgt_base_rects, x_center, vas_tgt_y_coords(2));
vas_tgt_rect_coords(:, 2) = CenterRectOnPointd(tgt_base_rects, x_center, y_center);

%%%%%% feedback info
feeds_vals = [0 50 500];
feed_col_high = [0 255 0];
feed_col_med = [255 191 0];
error_col = [255 0 0];
show_points = [1:5 11:15 21:25 31:35 41:45];
decay_rate = -4;
high_points = 5000;
low_points = 100;

%%%%%%% define fonts
Screen('TextFont', w, 'Courier New');
Screen('TextStyle', w, 1);
Screen('TextSize', w, 30);

%%%%%%%%%% define timings
time.before_first_trial = 2;
time.ramp = 1*expand_time;
time.fix_time = 1*expand_time;
time.pre_feed = .25*expand_time;
time.feedback = .5*expand_time;
time.cue_dur = .3*expand_time;
time.pre_cue_min = .4*expand_time;
time.pre_cue_max = .5*expand_time;
time.rest = 1*expand_time;
time.rest_total = 5;
time.tgt_on = .100*expand_time;
time.cue_tgt_soa = .05;
%time.dst_on = .100*expand_time;
time.resp_period = 2.5*expand_time;
%%%%%%%%%% define frames
frames.before_first_trial = round(time.before_first_trial/ifi);
frames.ramp = round(time.ramp/ifi);
frames.pre_feed = round(time.pre_feed/ifi);
frames.feedback = round(time.feedback/ifi);
frames.cue_dur = round(time.cue_dur/ifi);
frames.rest = round(time.rest/ifi);
frames.tgt_on = round(time.tgt_on/ifi);
frames.cue_tgt_soa = round(time.cue_tgt_soa/ifi);
frames.resp_period = round(time.resp_period/ifi);
frames.wait_frames = 1;

%%%%%%% define responses
if sess.response_order == 1
    task.responses = [KbName('v') KbName('g')];
elseif sess.response_order == 2
    task.responses = [KbName('g') KbName('v')];
end
move.vas = KbName('space');

%%%%%%%%%%%%% Initialise Eyelink
if sess.eye_on
    el=EyelinkInitDefaults(w); %Initialization
    if ~EyelinkInit(0) %Fail nice if you are going to fail
        fprintf('Eyelink Init aborted.\n');
        
        return;
    end
    [sess.v sess.vs]=Eyelink('GetTrackerVersion'); % optional
    % fprintf('Running experiment on a ''%s'' tracker.\n', vs ); % optional
    edfFile=sprintf('KGs%dv4.edf', sess.sub_num);
    Eyelink('Openfile', edfFile); %Create and open your Eyelink File
    % set calibration type.
    Eyelink('command', 'calibration_type = HV9'); %This is a typical 9 point calibration
    Eyelink('command', 'saccade_velocity_threshold = 35'); %Default from Eyelink Demo
    Eyelink('command', 'saccade_acceleration_threshold = 9500'); %Default from Eyelink Demo
    % set EDF file contents
    Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON'); %Event data to collect
    Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS'); %Sample data to collect
    el.backgroundcolour = screen_grey;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN EXPERIMENT
HideCursor;
count_trials = 0;
t_points = 0;
count_blocks = 0;
grand_trial_count = 0;

task.trial_structures = 5; % n infos for each trial 
task.t_trials_per_block = 200;  % 50 trials for orienting w/out value, 150 trials for value assoc + orient
task.trials = zeros(task.trial_structures, task.t_trials_per_block, length(sess.blocks));
task.trials_per_block = 50; %%%%% this is actually break (not block)
task.all_trials = task.t_trials_per_block*length(sess.blocks);
task.t_blocks = task.all_trials/task.trials_per_block;

%%%%%%%%%%%%% trial info to collect
trials.cor_resp = zeros(1,task.t_trials_per_block,length(sess.blocks));
trials.resp = zeros(1,task.t_trials_per_block,length(sess.blocks));
trials.rts = zeros(1,task.t_trials_per_block,length(sess.blocks));
trials.stim_on = zeros(1,task.t_trials_per_block,length(sess.blocks));
trials.cue_on = zeros(1,task.t_trials_per_block,length(sess.blocks));
trials.soa_on = zeros(1,task.t_trials_per_block,length(sess.blocks));
trials.target_on = zeros(1,task.t_trials_per_block,length(sess.blocks));
trials.rew_rec = zeros(1,task.t_trials_per_block,length(sess.blocks));
trials.distances = zeros(1, task.t_trials_per_block/task.trials_per_block, length(sess.blocks));
time.pre_cue_this_trial = zeros(1,task.t_trials_per_block,length(sess.blocks));

% starting instructions
Screen('TextStyle', w, 1);
Screen('TextSize', w, 30);
instructions = 'press any key to start';
DrawFormattedText(w, instructions, 'Center', 'Center', black, 115);
Screen('Flip', w);
while (KbCheck); end; while (~KbCheck); end

count_rewards = 0;

too_fast = 0;
for prob_types = 1:length(sess.blocks)
    
    %%%%%%%%%%%%% trial structure
    %%%%%%%%%% vals = 1 = high, 2 = low
    %%%%%%%%%% dirs = 1 = left, 2 = right
    %%%%%%%%%% reward = 1 = reward if correct, 0 = no reward
    %%%%%%%%%% left val, right val, cue dir, tgt loc, reward
    
    % get fixation
    alert_fix_coords = fix_all{1};
    alert_base_cols = alert_base_cols_all{1};
    alert_fix_cols = alert_fix_cols_all{1};
    
    trial_types = get_trials_for_block(  sess.blocks( prob_types ) );
    task.trials(:,:,prob_types) = trial_types;  % save the trials
    t_trials = task.t_trials_per_block; % set total number of trials for that block
    
    count_vas = 0;
    
    static = mod( sess.blocks ( prob_types ), 2 );
    
    for count_trials = 1:t_trials
   
    grand_trial_count = grand_trial_count + 1;
    
        if grand_trial_count == 1
            if sess.eye_on
                EyelinkDoTrackerSetup(el);
                %Manually do initial calibration, will come back to code here
                WaitSecs(0.5);
                
            end
        end
         
        if count_trials == 1
            
            % instructions for no points
            Screen('TextStyle', w, 1);
            Screen('TextSize', w, 30);
            instructions = 'no points this time...\nbut still be as accurate and as quick\n as you can!\npress any key to continue';
            DrawFormattedText(w, instructions, 'Center', 'Center', black, 115);
            Screen('Flip', w);
            while (KbCheck); end; while (~KbCheck); end
            current_cols = stooge_col;
            
        elseif count_trials == 41 
            
            if static == 1
                instructions = 'this time the points stays the same all the time\nstill be as accurate and as quick\n as you can!\npress any key to continue\ngood luck!';
            elseif static == 0
                instructions = 'this time the points run out - time counts!\nstill be as accurate and as quick\n as you can!\npress any key to continue\ngood luck!';
            end
            % instructions for points
            Screen('TextStyle', w, 1);
            Screen('TextSize', w, 30);
            
            DrawFormattedText(w, instructions, 'Center', 'Center', black, 115);
            Screen('Flip', w);
            while (KbCheck); end; while (~KbCheck); end
            current_cols = task.shape_col;
            
        end
        
        %%%%%% get each trial stimulus
        left_col = current_cols(trial_types(1,count_trials),:);
        right_col = current_cols(trial_types(2,count_trials),:);
        cue_dir = trial_types(3,count_trials);
        cue_cols = alert_fix_cols{trial_types(3,count_trials)};
        tgt_loc = trial_types(4,count_trials);
        if tgt_loc == 1
            tgt_x = target_center_left-x_offset;
            dst_x = target_center_right-x_offset;
            
        else
            tgt_x = target_center_right-x_offset;
            dst_x = target_center_left-x_offset;
        end
        c_tgt_id = randperm(length(tgts),1);
        c_tgt = tgts{c_tgt_id};
        trials.cor_resp(1,count_trials,prob_types) = c_tgt_id;
        c_dst_id =  randperm(length(dsts),1);
        c_dst = dsts{c_dst_id};
        
        %%%%%%%%%%% timing
        time.pre_cue_this_trial(1,count_trials,prob_types) = (time.pre_cue_max - time.pre_cue_min)*rand(1)+time.pre_cue_min;
        fix_frames_this_trial = round(time.pre_cue_this_trial(1,count_trials,prob_types)/ifi);
        cue_frames_this_cycle = frames.cue_dur;
        
        reward_trial = trial_types(5,count_trials); %%%%% correct
        
        if sess.eye_on
            % Start the eyetracker before each trial
            Eyelink('StartRecording'); % Start recording
            Eyelink('Message', 'SYNCTIME'); %Send a message to mark trial start
            trialStartTimeE(grand_trial_count) = Eyelink('TrackerTime'); %The time since you started the tracker (Optional)
            trialStartTimeOffset(grand_trial_count) = Eyelink('TimeOffset'); %The difference between PC and tracker time (Optional)
            Eyelink('Message', 'TRIALID %d', grand_trial_count); % Send a message with the trial number
            
            %send interest areas online (you can do this offline to!)
            Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, tgt_rect_coords(1,1),tgt_rect_coords(2,1),tgt_rect_coords(3,1),tgt_rect_coords(4,1), 'left'); %Set up the interest areas (Optional)
            Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 2, tgt_rect_coords(1,2),tgt_rect_coords(2,2),tgt_rect_coords(3,2),tgt_rect_coords(4,2), 'right');%Set up the interest areas (Optional)
            Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 3, alert_fix_check_coords(1), alert_fix_check_coords(2), alert_fix_check_coords(3), alert_fix_check_coords(4), 'fix');
            
        end
        
        if sess.eye_on
            
            Eyelink('Message', 'StimOnset');
        end
        
        %%%%%%%% present fix and stimuli
        Priority(topPriorityLevel);
        vbl = Screen('Flip', w);
        trials.stim_on(1,count_trials,prob_types) = GetSecs;
        for frame = 1:fix_frames_this_trial

            Screen( 'DrawLines', w, alert_fix_coords, alert_fix_width, alert_base_cols, [x_center, y_center], 1 );
            % calling fill oval twice as these will be flickered in later versions
            Screen('FrameOval', w, left_col, tgt_rect_coords(:, 1), oval_width); %%%%% dd - this should be green and left
            Screen('FrameOval', w, right_col, tgt_rect_coords(:, 2), oval_width);
            Screen('DrawingFinished', w);
            vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
            
        end
        
        if sess.eye_on
            
            Eyelink('Message', 'CueOnset');
        end
        vbl = Screen('Flip', w);
        trials.cue_on(1,count_trials,prob_types) = GetSecs;
        for frame = 1:cue_frames_this_cycle
            
            Screen( 'DrawLines', w, alert_fix_coords, alert_fix_width, cue_cols, [x_center, y_center], 1 );
            Screen('FrameOval', w, left_col, tgt_rect_coords(:, 1), oval_width); %%%%% dd - this should be green and left
            Screen('FrameOval', w, right_col, tgt_rect_coords(:, 2), oval_width);
            Screen('DrawingFinished', w);
            vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
        end
        
        eye_test = 0;    % check eyes in the right place
        if sess.eye_on
            if Eyelink('NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                evt = Eyelink('NewestFloatSample');
                
                % if we do, get current gaze position from sample
                x = evt.gx(2); % +1 as we're accessing MATLAB array
                y = evt.gy(2);
                % do we have valid data and is the pupil visible?
                if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(2)>0
                    eyemx=x;
                    eyemy=y;
                    erroreye=0;
                    % or are the eyes outside of fixation
                    if  (eyemx< alert_fix_check_coords(1) || eyemx > alert_fix_check_coords(3) || eyemy< alert_fix_check_coords(2) || eyemy> alert_fix_check_coords(4))
                        
                        for frame = 1:frames.rest
                            DrawFormattedText(w, 'Too Fast!', 'Center', 'Center', [255 0 0], 115);
                            vbl = Screen('Flip', w);
                        end
                        too_fast = too_fast + 1;
                        %trials.rem_trials(:,too_fast) = task.trial_struct(:,count_trials); % if this works, will create a script that runs only these trials
                        %%%% WORK OUT HOW TO LOG FILES/MAKE PEOPLE DO MORE
                        eye_test = 1;
                        
                    end
                end
            end
        end
        
        if eye_test
            trials.target_on( 1,count_trials,prob_types ) = NaN;
            trials.resp( 1,count_trials,prob_types ) = NaN;
            trials.rts( 1,count_trials,prob_types ) = NaN;
        else
            %%%%%%% present cue-tgt soa
            vbl = Screen('Flip', w);
            trials.soa_on(1,count_trials,prob_types) = GetSecs;
            for frame = 1:frames.cue_tgt_soa
                
            Screen( 'DrawLines', w, alert_fix_coords, alert_fix_width, alert_base_cols, [x_center, y_center], 1 );
                % calling fill oval twice as these will be flickered in later versions
                Screen('FrameOval', w, left_col, tgt_rect_coords(:, 1), oval_width); %%%%% dd - this should be green and left
                Screen('FrameOval', w, right_col, tgt_rect_coords(:, 2), oval_width);
                Screen('DrawingFinished', w);
                vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
            end
            
            if sess.eye_on
                
                Eyelink('Message', 'TgtOn');
            end
            
            %%%%%%%% present target for 100 ms
            Screen('TextStyle', w, 0);
            Screen('TextSize', w, 20);
            %         Priority(topPriorityLevel);
            vbl = Screen('Flip',w);
            key_down = 0;
            trials.target_on(1,count_trials,prob_types) = GetSecs;
            for frame = 1:frames.tgt_on
                
                Screen( 'DrawLines', w, alert_fix_coords, alert_fix_width, alert_base_cols, [x_center, y_center], 1 );
                Screen('FrameOval', w, left_col, tgt_rect_coords(:, 1), oval_width); %%%%% dd - this should be green and left
                Screen('FrameOval', w, right_col, tgt_rect_coords(:, 2), oval_width);
                Screen('DrawText', w, c_tgt, tgt_x, tgt_y, target_col);
                Screen('DrawText', w, c_dst, dst_x, dst_y, target_col);
                Screen('DrawingFinished', w);
                vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
                
                if dbug% if dbug mode do some dummy responses
                    trials.resp(1,count_trials,prob_types) = trials.cor_resp(1,count_trials,prob_types);%trials.cor_resp(1,count_trials,prob_types);
                    trials.rts(1,count_trials,prob_types) = rand(1);
                    key_down = 1;
                else
                    
                    %%%%%% add response collection code here (polling for quick
                    %%%%%% responses here
                    if ~key_down
                        [key_is_down,secs,key_code] = KbCheck;
                        if any(key_code(task.responses))
                            rt_stamp = GetSecs;
                            if key_code(task.responses(1))
                                trials.resp(1,count_trials,prob_types) = 1;
                            elseif key_code(task.responses(2))
                                trials.resp(1,count_trials,prob_types) = 2;
                            end
                            trials.rts(1,count_trials,prob_types) = rt_stamp - trials.target_on(1,count_trials,prob_types);
                            key_down = 1;
                        end
                    end
                end
                
            end
            
            % collect response
            vbl = Screen('Flip', w);
            for frame = 1:frames.resp_period
                
                Screen( 'DrawLines', w, alert_fix_coords, alert_fix_width, alert_base_cols, [x_center, y_center], 1 );
                Screen('FrameOval', w, left_col, tgt_rect_coords(:, 1), oval_width);
                Screen('FrameOval', w, right_col, tgt_rect_coords(:, 2), oval_width);
                Screen('DrawingFinished', w);
                vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
                
                %%%%%% add response collection code here (polling for quick
                %%%%%% responses here
                if ~key_down
                    [key_is_down,secs,key_code] = KbCheck;
                    if any(key_code(task.responses))
                        rt_stamp = GetSecs;
                        if key_code(task.responses(1))
                            trials.resp(1,count_trials,prob_types) = 1;
                        elseif key_code(task.responses(2))
                            trials.resp(1,count_trials,prob_types) = 2;
                        end
                        trials.rts(1,count_trials,prob_types) = rt_stamp - trials.target_on(1,count_trials,prob_types);
                        key_down = 1;
                        break
                    end
                end
            end
            
            if sess.eye_on
                
                Eyelink('Message', 'PreFeedOn');
            end
            %%%%%%%% present pre-reward stim
            vbl = Screen('Flip', w);
            for frame = 1:frames.pre_feed
                
                Screen( 'DrawLines', w, alert_fix_coords, alert_fix_width, alert_base_cols, [x_center, y_center], 1 );
                Screen( 'FrameOval', w, left_col, tgt_rect_coords(:, 1), oval_width); %%%%% dd - this should be green and left
                Screen( 'FrameOval', w, right_col, tgt_rect_coords(:, 2), oval_width);
                Screen( 'DrawingFinished', w);
                vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
                
            end
            
            %%%%%% if still no response assign NaN
            if ~key_down
                trials.resp(1,count_trials,prob_types) = NaN;
                trials.rts(1,count_trials,prob_types) = NaN;
            end
            
            %%%%%% is the response incorrect?
            if trials.cor_resp(1,count_trials,prob_types) ~= trials.resp(1,count_trials,prob_types)
                
                feed_case = 2;
                feed_col_this_trial = error_col;
                reward_fix_width = 4;
                
            else %%%%%% for correct trials
                
                if any(reward_trial)
                    feed_case = 1;
                    % get the value of the target location
                    % a) colour locations
                    tmp_locs = [trial_types(1,count_trials) trial_types(2,count_trials)];
                    % a, did target appear right or left
                    tmp_reward = tmp_locs(tgt_loc);
                    
                    if tmp_reward == 1
                        
                        if static == 1
                            trials.rew_rec(1,count_trials,prob_types) = high_points;
                        elseif static == 0
                            trials.rew_rec(1,count_trials,prob_types) = round(high_points*exp(decay_rate*trials.rts(1,count_trials,prob_types)));
                        end
                        feed_col_this_trial = feed_col_high;
                        
                        reward_text = num2str(trials.rew_rec(1,count_trials,prob_types));
                        count_rewards = count_rewards + trials.rew_rec(1,count_trials,prob_types);
                        
                    elseif tmp_reward == 2
                        
                        if static == 1
                            trials.rew_rec(1,count_trials,prob_types) = low_points;
                        elseif static == 0
                            trials.rew_rec(1,count_trials,prob_types) = round(low_points*exp(decay_rate*trials.rts(1,count_trials,prob_types)));
                        end
                        feed_col_this_trial = feed_col_med;
                        reward_text = num2str(trials.rew_rec(1,count_trials,prob_types));
                        count_rewards = count_rewards + trials.rew_rec(1,count_trials,prob_types);
                    end
                else
                    
                    feed_case = 2;
                    trials.rew_rec(1,count_trials,prob_types) = NaN; %%%%%%%%%% or no point available and standard fix during feedback
                    feed_col_this_trial = alert_base_cols;
                    reward_fix_width = 4;
                end
            end
            
            Screen('TextStyle', w, 1);
            Screen('TextSize', w, 20);
            
            vbl = Screen('Flip', w);
            for frame = 1:frames.feedback
                
                if feed_case == 2
                    Screen('DrawLines', w, alert_fix_coords, reward_fix_width, feed_col_this_trial, [x_center, y_center], 1);

                else
                    DrawFormattedText(w, reward_text, 'Center', 'Center', feed_col_this_trial, 115);
                    
                end
                Screen('FrameOval', w, neut_col, tgt_rect_coords(:,1), oval_width/2);
                Screen('FrameOval', w, neut_col, tgt_rect_coords(:,2), oval_width/2);
                Screen('DrawingFinished', w);
                vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
            end
        end
        Priority(0);
        
        %%%%%%%%% save stuff - wasn't working on PC.
        if count_trials == 1 && prob_types == 1
            save(save_mat, 'frames', 'time', 'sess', 'task', 'trials');
        else
            save(save_mat, 'task', 'trials', '-append');
        end
        %%% every N trials (as defined above), stop for a break
        
        if ~mod( grand_trial_count, task.trials_per_block)
            
            % get the vas measure
            vbl = Screen('Flip', w);
            SetMouse( x_center, y_center, w );
            offsetSet = 0;
            space_press = 0;
            move_rect = vas_tgt_rect_coords(:, 1);
            count_vas = count_vas + 1;
            sx = x_center;
            sy = y_center - 150;
            while ~space_press
                
                DrawFormattedText(w, 'put the target on the line\npress space when done', x_center - 200, y_center - 250, [0 0 255], 115);
                Screen('DrawLines', w, vas_coords, vas_width_pix, white, [x_center, y_center], 1);
                Screen('FrameOval', w, left_col, vas_tgt_rect_coords(:, 2), oval_width);
                
                
                [mx, my, buttons] = GetMouse( w ); % get mouse coords
                [cx, cy] = RectCenter(vas_tgt_rect_coords(:,1)); % get central position of moving oval
                inside = IsInRect(mx, my, move_rect); % is the mouse inside the target?
                
                if inside == 1 && sum(buttons) > 0 && offsetSet == 0
                    % If the mouse cursor is inside the square and a mouse button is being
                    % pressed and the offset has not been set, set the offset and signal
                    % that it has been set
                    dx = mx - cx;
                    dy = my - cy;
                    offsetSet = 1;
                end
                
                % If we are clicking on the square allow its position to be modified by
                % moving the mouse, correcting for the offset between the centre of the
                % square and the mouse position
                if inside == 1 && sum(buttons) > 0
                    sx = mx - dx;
                    sy = my - dy;
                end
                
                % centre oval on new position
                move_rect = CenterRectOnPointd(tgt_base_rects, sx, sy);
                Screen('FrameOval', w, right_col, move_rect, oval_width); % draw oval
                % Draw a white dot where the mouse cursor is
                Screen('DrawDots', w, [mx, my], 10, white, [], 2);
                % Flip to the screen
                vbl  = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
                % Check to see if the mouse button has been released and if so reset
                % the offset cue
                if sum(buttons) <= 0
                    offsetSet = 0;
                end
                % check for space bar
                if ~space_press
                    [key_is_down,secs,key_code] = KbCheck;
                    if key_code( move.vas )
                        
                        trials.distances(1, count_vas, prob_types) = sx - x_center;
                        space_press = 1;
                        WaitSecs(1);
                    end
                end
            end
            
            
            if count_trials == t_trials && prob_types == length(sess.blocks)
            else
                t_points = count_rewards;
                Screen('TextStyle', w, 1);
                Screen('TextSize', w, 30);
                count_blocks = count_blocks + 1;
                points_info = sprintf( 'you have earned %d points!!', count_rewards );
                instructions = sprintf('woo-hoo\n! %d/%d blocks complete\n\n take a break \npress any key when you are ready', count_blocks, task.t_blocks);
                DrawFormattedText(w, points_info, x_center - 175, y_center - 250, [0 255 0], 115);
                DrawFormattedText(w, instructions, 'Center', 'Center', black, 115);
                Screen('Flip', w);
                
                if dbug
                else
                    while (KbCheck); end; while (~KbCheck); end
                end
                
                if dbug
                else
                    while (KbCheck); end; while (~KbCheck); end
                end
                if sess.eye_on
                    EyelinkDoTrackerSetup(el);
                    %Manually do initial calibration, will come back to code here
                    WaitSecs(0.5);
                    
                    
                end
                
                
                %%%%% preparation screen %%%%%%%%%%%%% check variables on this
                instructions = 'next trial starting in...\n';
                vbl = Screen('Flip', w);
                for i = 1:time.rest_total
                    
                    for frame = 1:frames.rest
                        instruct = [instructions num2str(time.rest_total+1 - i)];
                        DrawFormattedText(w, instruct, 'Center', 'Center', black, 115);
                        Screen('Flip', w);
                    end
                end
            end
        end
        
    end
end


if dbug
else
    while (KbCheck); end; while (~KbCheck); end
end

if sess.eye_on == 1
    Eyelink( 'StopRecording' );
    Eyelink( 'CloseFile' );
    Eyelink( 'ReceiveFile', upper( edfFile ));
end


t_points = count_rewards;
final_text = 'Hooray - you are finished!';
Screen('TextStyle', w, 1);
Screen('TextSize', w, 30);
points_info = sprintf('you have earned %d points!', t_points);
DrawFormattedText(w, points_info, x_center - 175, y_center - 250, [0 255 0], 115);
DrawFormattedText(w, final_text, 'Center', 'Center', black, 115);
vbl = Screen('Flip', w);
while (KbCheck); end; while (~KbCheck); end
Screen('CloseAll');

