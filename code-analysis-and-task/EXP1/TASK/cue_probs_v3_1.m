%%%%%%%%%%%%% K. Garner - Sept 2016, University of Birmingham
%%%%%%%%%%%%% task is a relative value cueing task where participants respond to a cued
%%%%%%%%%%%%% target (.98/.02,.95/.05, .9/.1, .8/.2, .6/.4), irrespective of
%%%%%%%%%%%%% the reward that target's placeholder represents. 
%%%%%%%%%%%%% Cue-target SOA = 50 ms (Lanthier et al, 2015).
%%%%%%%%%%%%% tgts on for 100 ms

% _______________________________________________________________________________________________________________
% DEPENDENCIES
% coded on Matlab 2013a using PTB v3.0.12
% run on Stone SOFREP-144
% with Asus VG278HE

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
sess.sub_num = input('Subject Number? ');
sess.session = input('Session? ');
sess.blocks = [1 2 3 4 5]; % .98/.02, .94/.06, .9/.1, .8/.2, .6/.4
sess.blocks = sess.blocks(randperm(length(sess.blocks)));
sess.response_order = input('Response Order? enter 1 or 2 ');
sess.shapes = input('Shape Order? enter 1 or 2 ');
sess.prev_points = input('Points? ');
%%%%%%%%% set time cheats
expand_time = 1;

% set randomisation seed based on sub/sess number
r_num = [num2str(sess.sub_num) num2str(sess.session)];
r_num = str2double(r_num);
rand('state',r_num);
randstate = rand('state');

save_mat = ['KG_exp3_cueprob_v3_1_' num2str(sess.sub_num) '_s' num2str(sess.session)];

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
[w, rect] = Screen('OpenWindow', screen_nums); % ,...
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
tmp = [86.7 74.562 80.631; 119.85 86.2321 1.1985]; % check coulour luminance
task.shape_col = [tmp(sess.shapes(1),:); tmp(sess.shapes(2),:)]; % now call 1 for no reward, 2 for high
oval_width = 8;
neut_col = [155 155 155];
%%%%%%%% now draw the targets, which can be an upright or inverted cross
target_col = [90 90 90];
target_center_left = x_center - center_offset; %%%%%% place to center target on
target_center_right = x_center + center_offset;
tgts = {'H','N'};
dsts = {'K','Z'};
msks = {'8','8'};
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
base_cols_a = [200 200 200; 200 200 200; 200 200 200; 200 200 200; 200 200 200; 200 200 200; 200 200 200; 200 200 200]';
base_cols_b = [200 200 200; 200 200 200; 200 200 200; 200 200 200]';
alert_base_cols_all = {base_cols_a base_cols_b};
alert_base_cols_all = {alert_base_cols_all{fix_allocate(1)} alert_base_cols_all{fix_allocate(2)}}; 
alert_fix_cols_a = {[50 50 50; 50 50 50; 50 50 50; 50 50 50; 200 200 200; 200 200 200; 200 200 200; 200 200 200]',[200 200 200; 200 200 200; 200 200 200; 200 200 200; 50 50 50; 50 50 50; 50 50 50; 50 50 50]'};
alert_fix_cols_b = {[50 50 50; 50 50 50; 200 200 200; 200 200 200]',[200 200 200; 200 200 200; 50 50 50; 50 50 50]'};
alert_fix_cols_all = {alert_fix_cols_a alert_fix_cols_b};
alert_fix_cols_all = {alert_fix_cols_all{fix_allocate(1)} alert_fix_cols_all{fix_allocate(2)}};
alert_fix_width = 4;
reward_fix_width = 8;
alert_fix_check_coords = CenterRectOnPointd([0 0 100 100], x_center, y_center);

%%%%%% feedback info
feeds_vals = [0 50 500];
feed_col_high = [0 255 0];
feed_col_med = [255 191 0];
error_col = [255 0 0];
show_points = [1:5 11:15 21:25 31:35 41:45];

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
% frames.fix_time = round(time.fix_time/ifi);
frames.pre_feed = round(time.pre_feed/ifi);
frames.feedback = round(time.feedback/ifi);
frames.cue_dur = round(time.cue_dur/ifi);
frames.rest = round(time.rest/ifi);
frames.tgt_on = round(time.tgt_on/ifi);
%frames.dst_on = round(time.dst_on/ifi);
frames.cue_tgt_soa = round(time.cue_tgt_soa/ifi);
frames.resp_period = round(time.resp_period/ifi);
frames.wait_frames = 1;


%%%%%%% define responses
if sess.response_order == 1
    task.responses = [KbName('v') KbName('g')];
elseif sess.response_order == 2
    task.responses = [KbName('g') KbName('v')];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN EXPERIMENT
HideCursor;
count_trials = 0;
t_points = 0;
count_blocks = 0;
grand_trial_count = 0;

task.trial_structures = 5;
task.t_trials_per_block = 100;
task.trials = zeros(task.trial_structures,task.t_trials_per_block,length(sess.blocks));
task.trials_per_block = 25; %%%%% this is actually break (not block)
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
%trials.distracter_on = zeros(1,task.t_trials_per_block,length(sess.blocks));
trials.rew_rec = zeros(1,task.t_trials_per_block,length(sess.blocks));
time.pre_cue_this_trial = zeros(1,task.t_trials_per_block,length(sess.blocks));
% starting instructions
Screen('TextStyle', w, 1);
Screen('TextSize', w, 30);
instructions = 'find the target\npress any key to start';
DrawFormattedText(w, instructions, 'Center', 'Center', black, 115);
Screen('Flip', w);
while (KbCheck); end; while (~KbCheck); end

count_rewards = sess.prev_points;


for prob_types = 1:length(sess.blocks)
    
    %%%%%%%%%%%%% trial structure
    %%%%%%%%%% vals = 1 = high, 2 = low
    %%%%%%%%%% dirs = 1 = left, 2 = right
    %%%%%%%%%% reward = 1 = reward if correct, 0 = no reward
    %%%%%%%%%% left val, right val, cue dir, tgt loc, reward

    alert_fix_coords = fix_all{1};
    alert_base_cols = alert_base_cols_all{1};
    alert_fix_cols = alert_fix_cols_all{1};
    
    if sess.blocks(prob_types) == 1
        cue_prob = [24, 1];
    elseif sess.blocks(prob_types) == 2
        cue_prob = [23, 2];
    elseif sess.blocks(prob_types) == 3
        cue_prob = [22, 3];
    elseif sess.blocks(prob_types) == 4
        cue_prob = [20, 5];
    elseif sess.blocks(prob_types) == 5
        cue_prob = [15, 10];
    end
    
    %%%%%%%%%%%%%%%% hvl conditions
    %%%%%% left cue
    hvl_leftcue_val = repmat([1, 2, 1, 1], [cue_prob(1),1]);
    hvl_leftcue_val(1:(round(cue_prob(1)*.8)),5) = 1; %%%%%% reward available on 80 % of trials
    hvl_leftcue_inv = repmat([1, 2, 1, 2], [cue_prob(2),1]);
    hvl_leftcue_inv(1:(round(cue_prob(2)*.8)),5) = 1;
    
    %%%%%% right cue
    hvl_rightcue_val = repmat([2, 1, 2, 2],[cue_prob(1),1]);
    hvl_rightcue_val(1:(round(cue_prob(1)*.8)),5) = 1; %%%%%% reward available on 80 % of trials
    hvl_rightcue_inv = repmat([2, 1, 2, 1],[cue_prob(2),1]);
    hvl_rightcue_inv(1:(round(cue_prob(2)*.8)),5) = 1;
    
    %%%%%%%%%%%%%%% lvh conditions
    %%%%%% left cue
    lvh_leftcue_val = repmat([2, 1, 1, 1], [cue_prob(1),1]);
    lvh_leftcue_val(1:(round(cue_prob(1)*.8)),5) = 1;
    lvh_leftcue_inv = repmat([2, 1, 1, 2], [cue_prob(2),1]);
    lvh_leftcue_inv(1:(round(cue_prob(2)*.8)),5) = 1; %%%%%% reward available on 80 % of trials
    
    %%%%%% right cue
    lvh_rightcue_val = repmat([1, 2, 2, 2], [cue_prob(1),1]);
    lvh_rightcue_val(1:(round(cue_prob(1)*.8)),5) = 1;
    lvh_rightcue_inv = repmat([1, 2, 2, 1], [cue_prob(2),1]);
    lvh_rightcue_inv(1:(round(cue_prob(2)*.8)),5) = 1; %%%%%% reward available on 80 % of trials
    
    trial_types = [hvl_leftcue_val; hvl_leftcue_inv; hvl_rightcue_val; hvl_rightcue_inv;...
        lvh_leftcue_val; lvh_leftcue_inv; lvh_rightcue_val; lvh_rightcue_inv]';
    t_trials = task.t_trials_per_block;
    trial_types = trial_types(:,randperm(t_trials));
    task.trials(:,:,prob_types) = trial_types;
    
    for count_trials = 1:t_trials
        
        grand_trial_count = grand_trial_count + 1;
        
        %%%%%% get each trial stimulus
        left_col = task.shape_col(trial_types(1,count_trials),:); %%%% during debug (dd)
        right_col = task.shape_col(trial_types(2,count_trials),:); %%%%%%% (dd)
        cue_dir = trial_types(3,count_trials); %%%%%
        cue_cols = alert_fix_cols{trial_types(3,count_trials)}; %%%%%
        tgt_loc = trial_types(4,count_trials); %%%%% where should the target go?
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
        
        %%%%%%%% present fix and stimuli
        Priority(topPriorityLevel);
        vbl = Screen('Flip', w);
        trials.stim_on(1,count_trials,prob_types) = GetSecs;
        for frame = 1:fix_frames_this_trial
            Screen('DrawLines', w, alert_fix_coords, alert_fix_width, alert_base_cols, [x_center, y_center], 1);
            % calling fill oval twice as these will be flickered in later versions
            Screen('FrameOval', w, left_col, tgt_rect_coords(:, 1), oval_width); %%%%% dd - this should be green and left
            Screen('FrameOval', w, right_col, tgt_rect_coords(:, 2), oval_width);
            Screen('DrawingFinished', w);
            vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
            
        end
%         Priority(0);
        
        %%%%%%%% present cue stimuli
%         Priority(topPriorityLevel);
        vbl = Screen('Flip', w);
        trials.cue_on(1,count_trials,prob_types) = GetSecs;
        for frame = 1:cue_frames_this_cycle
            
            Screen('DrawLines', w, alert_fix_coords, alert_fix_width, cue_cols, [x_center, y_center], 1);
            Screen('FrameOval', w, left_col, tgt_rect_coords(:, 1), oval_width); %%%%% dd - this should be green and left
            Screen('FrameOval', w, right_col, tgt_rect_coords(:, 2), oval_width);
            Screen('DrawingFinished', w);
            vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
        end
%         Priority(0);

        %%%%%%% present cue-tgt soa
        vbl = Screen('Flip', w);
        trials.soa_on(1,count_trials,prob_types) = GetSecs;
        for frame = 1:frames.cue_tgt_soa
            
            Screen('DrawLines', w, alert_fix_coords, alert_fix_width, alert_base_cols, [x_center, y_center], 1);
            % calling fill oval twice as these will be flickered in later versions
            Screen('FrameOval', w, left_col, tgt_rect_coords(:, 1), oval_width); %%%%% dd - this should be green and left
            Screen('FrameOval', w, right_col, tgt_rect_coords(:, 2), oval_width);
            Screen('DrawingFinished', w);
            vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
        end
        
        
        %%%%%%%% present target for 100 ms
        Screen('TextStyle', w, 0);
        Screen('TextSize', w, 20);
%         Priority(topPriorityLevel);
        vbl = Screen('Flip',w);
        key_down = 0;
        trials.target_on(1,count_trials,prob_types) = GetSecs;
        for frame = 1:frames.tgt_on
            
            Screen('DrawLines', w, alert_fix_coords, alert_fix_width, alert_base_cols, [x_center, y_center], 1);
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
        
        %%%%%%%%% now present just the reward stim and poll for response
%         Priority(topPriorityLevel);
        vbl = Screen('Flip', w);
        for frame = 1:frames.resp_period
            Screen('DrawLines', w, alert_fix_coords, alert_fix_width, alert_base_cols, [x_center, y_center], 1);
            Screen('FrameOval', w, left_col, tgt_rect_coords(:, 1), oval_width); %%%%% dd - this should be green and left
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
        
        %%%%%%%% present pre-reward stim
%         Priority(topPriorityLevel);
        vbl = Screen('Flip', w);
        for frame = 1:frames.pre_feed
            Screen('DrawLines', w, alert_fix_coords, alert_fix_width, alert_base_cols, [x_center, y_center], 1);
            Screen('FrameOval', w, left_col, tgt_rect_coords(:, 1), oval_width); %%%%% dd - this should be green and left
            Screen('FrameOval', w, right_col, tgt_rect_coords(:, 2), oval_width);
            Screen('DrawingFinished', w);
            vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
            
        end
%         Priority(0);
        
        %%%%%% if still no response assign NaN
        if ~key_down
            trials.resp(1,count_trials,prob_types) = NaN;
            trials.rts(1,count_trials,prob_types) = NaN;
        end
        
        %%%%%% is the response incorrect?
        if trials.cor_resp(1,count_trials,prob_types) ~= trials.resp(1,count_trials,prob_types)
            feed_case = 2;
%             if any(reward_trial) % then scored the rewards received with a 0
%                 trials.rew_rec(1,count_trials,prob_types) = 0; % if error then get error feedback
%             else trials.rew_rec(1,count_trials,prob_types) = NaN; % if not score as nan because no reward was available
%             end
            feed_col_this_trial = error_col;
            reward_fix_width = 4;
        else %%%%%% for correct trials

                feed_case = 2;
                trials.rew_rec(1,count_trials,prob_types) = NaN; %%%%%%%%%% or no point available and standard fix during feedback
                feed_col_this_trial = alert_base_cols;
                reward_fix_width = 4;
        end
        
        Screen('TextStyle', w, 1);
        Screen('TextSize', w, 20);
%         Priority(topPriorityLevel);
        vbl = Screen('Flip', w);
        for frame = 1:frames.feedback
%             if feed_case == 1                
%                 DrawFormattedText(w, reward_text, 'Center', y_center-15, feed_col_this_trial, 115);
% 
%             else                
                Screen('DrawLines', w, alert_fix_coords, reward_fix_width, feed_col_this_trial, [x_center, y_center], 1);

%             end
            Screen('FrameOval', w, neut_col, tgt_rect_coords(:,1), oval_width/2);
            Screen('FrameOval', w, neut_col, tgt_rect_coords(:,2), oval_width/2);
            Screen('DrawingFinished', w);
            vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
        end
%         Priority(0);
        
        %%%%%% go back to grey circles before the next trial
%         Priority(topPriorityLevel);
        vbl = Screen('Flip', w);
        for frame = 1:frames.pre_feed
            Screen('DrawLines', w, alert_fix_coords, alert_fix_width, alert_base_cols, [x_center, y_center], 1);
            Screen('FrameOval', w, neut_col, tgt_rect_coords(:, 1), oval_width/2); %%%%% dd - this should be green and left
            Screen('FrameOval', w, neut_col, tgt_rect_coords(:, 2), oval_width/2);
            Screen('DrawingFinished', w);
            vbl = Screen('Flip', w, vbl + (frames.wait_frames - 0.5) * ifi);
            
        end
        Priority(0);        
        
        %%%%%%%%% save stuff - wasn't working on PC.
        if count_trials == 1 && prob_types == 1
            save(save_mat, 'frames', 'time', 'sess', 'task', 'trials');
        else
            save(save_mat, 'task', 'trials', '-append');
        end
        %%% every N trials (as defined above), stop for a break
        
        if ~mod(grand_trial_count,task.trials_per_block)
            if count_trials == t_trials && prob_types == length(sess.blocks)
            else
                t_points = count_rewards;
                Screen('TextStyle', w, 1);
                Screen('TextSize', w, 30);
                count_blocks = count_blocks + 1;
                %points_info = sprintf('you have earned %d points!!', t_points);
                instructions = sprintf('woo-hoo\n! %d/%d blocks complete\n\n take a break \npress any key when you are ready',count_blocks,task.t_blocks);
                %DrawFormattedText(w, points_info, x_center - 175, y_center - 250, [0 255 0], 115);
                DrawFormattedText(w, instructions, 'Center', 'Center', black, 115);
                Screen('Flip', w);
                
                if dbug
                else
                    while (KbCheck); end; while (~KbCheck); end
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


t_points = count_rewards;
final_text = 'Hooray - you are finished!';
Screen('TextStyle', w, 1);
Screen('TextSize', w, 30);
%points_info = sprintf('you have earned %d points!', t_points);
%DrawFormattedText(w, points_info, x_center - 175, y_center - 250, [0 255 0], 115);
DrawFormattedText(w, final_text, 'Center', 'Center', black, 115);
vbl = Screen('Flip', w);
while (KbCheck); end; while (~KbCheck); end
Screen('CloseAll');

