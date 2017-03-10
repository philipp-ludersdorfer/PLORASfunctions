function [tasks, conditions] = ploras_fmri_conditions(index,complex)
% PLORAS_FMRI_CONDITIONS(INDEX,COMPLEX) can be used to get the names of the 
% tasks and conditions of the PLORAS 3 fMRI paradigm. 
%
% When used without input arguments (i.e. ‘ploras_fmri_conditions’), the 
% function will display all task names at the Matlab command window. 
%
% INDEX (optional): a numerical vector indicating the indices of the tasks 
% to be displayed (e.g. 1 or [1,5] or 1:5).
%
% You can also save the task names by simply providing output variables:
% tasks = ploras_fmri_conditions
%
% COMPLEX (optional): if set to 1, the function will additionally save
% the condition names of the (selected) tasks. Usage:
% [tasks conditions] = ploras_fmri_conditions(1:5,1)
%
% Philipp Ludersdorfer (last modified 29/11/2016)

alltasks = {
    '1_Pic_Sem_Ass'
    '2_Pic_Name_2_Ob'
    '3_Pic_Name_Verb'
    '4_Pic_Name_Sent'
    '5_Aud_Sem_Ass'
    '6_Read_1_Ob'
    '7_Aud_Rep_1_Ob'
    '8_Pic_Name_1_Ob'
    '9_Colour_Name'
    '10_Aud_Sounds'
    '11_Read_Pseudo'
    '12_Aud_Rep_Pseudo'
    '13_Gender_Name'
    '14_Aud_Spell_1_Ob'
    };

allconditions = {
    {
    '1_SR_2_syll'
    '1_UR_2_syll'
    '1_SR_3+syll'
    '1_UR_3+syll'
    }
    {
    '2_new_set_2_syll'
    '2_new_set_3+syll'
    '2_repeat_set_2_syll'
    '2_repeat_set_3+syll'
    }
    {
    '3_new_set'
    '3_new_set'
    '3_repeat_set'
    '3_repeat_set'
    }
    {
    '4_new_set'
    '4_new_set'
    '4_repeat_set'
    '4_repeat_set'
    }
    {
    '5_SR_new_set'
    '5_UR_new_set'
    '5_SR_repeat_set'
    '5_UR_repeat_set'
    }
    {
    '6_new_set_1_syll'
    '6_new_set_2+syll'
    '6_repeat_set_1_syll'
    '6_repeat_set_2+syll'
    }
    {
    '7_new_set_1_syll'
    '7_new_set_2+syll'
    '7_repeat_set_1_syll'
    '7_repeat_set_2+syll'
    }
    {
    '8_new_set_1_syll'
    '8_new_set_2+syll'
    '8_repeat_set_1_syll'
    '8_repeat_set_2+syll'
    }
    {
    '9_block_1'
    '9_block_2'
    '9_block_3'
    '9_block_4'
    }
    {
    '10_novel_sounds_1st_pres'
    '10_novel_sounds_2nd_pres'
    '10_repeat_sounds_1st_pres'
    '10_repeat_sounds_2nd_pres'
    }
    {
    '11_new_set_1_syll'
    '11_new_set_2_syll'
    '11_repeat_set_1_syll'
    '11_repeat_set_2_syll'
    }
    {
    '12_new_set_1_syll'
    '12_new_set_2_syll'
    '12_repeat_set_1_syll'
    '12_repeat_set_2_syll'
    }
    {
    '13_short_male'
    '13_short_female'
    '13_long_male'
    '13_long_female'
    }
    {
    '14_early_consistent'
    '14_early_inconsistent'
    '14_late_consistent'
    '14_late_inconsistent'}
    };

if nargin < 1
    tasks = alltasks; 
elseif nargin > 0
    tasks = alltasks(index);
    if (nargin > 1) && (complex == 1)
        conditions = allconditions(index);
    end
end