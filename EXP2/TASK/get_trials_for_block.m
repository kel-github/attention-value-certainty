function [trial_types] = get_trials_for_block( prob_cond )
    % K. Garner, UoB, 2016

    if prob_cond == 1
        cue_prob = [24, 16];
    elseif prob_cond == 2
        cue_prob = [24, 16];
    elseif prob_cond == 3
        cue_prob = [32, 8];
    elseif prob_cond == 4
        cue_prob = [32, 8];
    end
    
    %%%%%%%%%%%%%%%% stooge conditions
    sto_leftcue_val = repmat([1, 2, 1, 1, 0], [cue_prob(1)/2,1]);
    sto_leftcue_inval = repmat([1, 2, 1, 2, 0], [cue_prob(2)/2,1]);   
    sto_rightcue_val = repmat([2, 1, 2, 2, 0], [cue_prob(1)/2,1]);
    sto_rightcue_inval = repmat([2, 1, 2, 1, 0],[cue_prob(2)/2,1]);
    
    stooge_trials = [ sto_leftcue_val; sto_leftcue_inval; sto_rightcue_val; sto_rightcue_inval ];
    %%%%%%%%%%%%%%% now shuffle the colour conditions
    [srows, ~] = size( stooge_trials );
    stooge_trials( :, [1 2] ) = stooge_trials( randperm(srows), [1 2]); 
    %%%%%%%%%%%%%% randomise all the trials
    stooge_trials = stooge_trials( randperm(srows), : )';
    
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
    [ ~ , t_trials ] = size( trial_types );
    trial_types = [stooge_trials trial_types(:,randperm(t_trials))];
    
end
    
    
    
    