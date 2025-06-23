function pipeline_03_finalize(study_info)

pipeline='NEARICA_behav_v3';

% Number of subjects
n_subjects=size(study_info.participant_info,1);

epoch_types={'eye','task'};

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    
    % Path containing subject data
    subject_dir=fullfile(study_info.data_dir, 'derivatives', pipeline, subj_id);
    subject_data_dir=fullfile(subject_dir, '04_rereferenced_data');
    
    for e=1:length(epoch_types)
        epoch_type=epoch_types{e};
        fname=sprintf('%s_rereferenced_data_%s.set',subj_id, epoch_type);
    
        if exist(fullfile(subject_data_dir,fname),'file')==2
        
            % Load data
            EEG=pop_loadset('filepath', subject_data_dir,...
                'filename', fname);   
        
            load(fullfile(subject_dir,'05_zapped_data',sprintf('zapline_%s.mat',epoch_type)));
            data=data.*1e6;
            EEG.data=permute(data,[2 3 1]);
            EEG = pop_editset(EEG, 'setname', sprintf('%s_%',subj_id, epoch_type));
            pop_saveset(EEG, 'filepath', fullfile(subject_dir, 'processed_data'),...
                'filename', sprintf('%s_%s.set',subj_id, epoch_type));
        
            fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
            saveas(fig, fullfile(subject_dir,sprintf('12-zapped_psd_%s.png', epoch_type)));
        end        
    end
    close all
end
