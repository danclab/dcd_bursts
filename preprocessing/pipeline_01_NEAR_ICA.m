%function pipeline_01_NEAR_ICA(study_info)

addpath('/home/bonaiuto/eeglab2021.1');
addpath('adjusted_adjust');

ext='.set';

% Channel location file
channel_locations = '/home/bonaiuto/dcd_ross/data/Standard-10-10-Cap33.ced';

% Initialize the filters
% High-pass frequency
highpass = 1;
% Low-pass frequency
lowpass  = 100;

% event/condition markers
epoch_labels_1={'EyesOpen', 'EyesClosed'};
event_markers_1 = {'E  3','E  5'};

% epoch length in seconds
epoch_length_1 = [0 30];

epoch_labels_2={'Kaleidoscope',	'ObserveFine', 'ObserveGross',	'ExecuteFine', 'ExecuteGross'};
event_markers_2 = {'E  9','E 17','E 18','E 33','E 34','E 40'};

% epoch length in seconds
epoch_length_2 = [0 10];

% lower and upper voltage threshold (in mV)
volt_threshold = [-150 150];

% list of frontal channels to check for epoch-level channel interpolation
% (see manuscript for detail)
frontal_channels = {'Fp1', 'Fp2'}; 

% Parameters for NEAR - Bad Channels Detection
% flat channels
isFlat  = 1; % flag variable to enable or disable Flat-lines detection method (default: 1)
flatWin = 5; % tolerance level in s(default: 5)

% LOF (density-based)
isLOF       = 1;  % flag variable to enable or disable LOF method (default: 1)
dist_metric = 'seuclidean'; % Distance metric to compute k-distance
thresh_lof  = 2.5; % Threshold cut-off for outlier detection on LOF scores
isAdapt = 10; % The threshold will be incremented by a factor of 1 if the given threshold detects more than xx %
%of total channels (eg., 10); if this variable left empty [], no adaptive thresholding is enabled.

% Periodogram (frequency based)
isPeriodogram = 0; % flag variable to enable or disable periodogram method (default: 0)
frange        = [1 20]; % Frequency Range in Hz
winsize       = 1; % window length in s
winov         = 0.66; % 66% overlap factor
pthresh       = 4.5; % Threshold Factor to predict outliers on the computed energy

% Parameters for NEAR- Bad Segments Correction/Rejection using ASR %
rej_cutoff = 13;   % A lower value implies severe removal (Recommended value range: 20 to 30)
rej_mode   = 'off'; % Set to 'off' for ASR Correction and 'on for ASR Removal (default: 'on')
add_reject = 'off'; % Set to 'on' for additional rejection of bad segments if any after ASR processing (default: 'off')

% Parameter for interpolation %
interp_type = 'spherical'; % other values can be 'v4'. Please refer to pop_interp.m for more details.

%% Initialize output variables
lof_flat_channels={};
lof_channels={};
lof_periodo_channels={};
% Bad channels identified using LOF
lof_bad_channels={};
% number of bad channel/s due to channel/s exceeding xx% of artifacted epochs
ica_preparation_bad_channels=[];
% length of data (in second) fed into ICA decomposition
length_ica_data=[];
% total independent components (ICs)
total_ICs=[];
% number of artifacted ICs
ICs_removed=[];
% number of epochs before artifact rejection
total_eye_epochs_before_artifact_rejection=[];
% number of epochs after artifact rejection
total_eye_epochs_after_artifact_rejection=[];
total_eye_channels_interpolated=[];
% number of epochs before artifact rejection
total_task_epochs_before_artifact_rejection=[];
% number of epochs after artifact rejection
total_task_epochs_after_artifact_rejection=[];
total_task_channels_interpolated=[];
asr_tot_samples_modified=[];
asr_change_in_RMS=[];

%% Loop over all data files
for s_idx=1:size(study_info.participant_info,1)
    % Get subject ID from study info
    subject=study_info.participant_info.participant_id{s_idx};
    
    % Where original raw data is located
    subject_raw_data_dir=fullfile(study_info.data_dir, 'raw');
    
    % Where to put processed (derived) data
    subject_output_data_dir=fullfile(study_info.data_dir, 'derivatives', 'NEARICA_behav_v3', subject);
    
    if exist([subject_output_data_dir filesep '01_filtered_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '01_filtered_data'])
    end

    if exist([subject_output_data_dir filesep '02_near_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '02_near_data'])
    end

    if exist([subject_output_data_dir filesep '03_ica_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '03_ica_data'])
    end

    if exist([subject_output_data_dir filesep '04_rereferenced_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '04_rereferenced_data'])
    end
    
    if exist([subject_output_data_dir filesep '05_zapped_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '05_zapped_data'])
    end
    
    if exist([subject_output_data_dir filesep 'processed_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep 'processed_data'])
    end
    
    fprintf('\n\n\n*** Processing subject %s ***\n\n\n', subject);
    
    %% Import data
%     data_file_name=sprintf('%s.vhdr',subject);
%     EEG = pop_loadbv(subject_raw_data_dir, data_file_name);
%     EEG = eeg_checkset(EEG);

    data_file_name = sprintf('%s.vhdr', subject);
    EEG = pop_loadbv(subject_raw_data_dir, data_file_name, [], []);
%     EEG.event = [];  % clear default events loaded from .vmrk
    EEG = eeg_checkset(EEG);
    
    % Import custom events from .cmrk file
    cmrk_filepath = fullfile(subject_raw_data_dir, sprintf('%s.cmrk', subject));
    EEG = import_cmrk_events(EEG, cmrk_filepath);

    EEG = eeg_checkset(EEG);
    
    %% Import channel locations
    EEG=pop_chanedit(EEG, 'lookup','/home/bonaiuto/eeglab2021.1/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
    EEG = eeg_checkset( EEG );
    
    % Check whether the channel locations were properly imported. The EEG
    % signals and channel numbers should be same.
    if size(EEG.data, 1) ~= length(EEG.chanlocs)
        error('The size of the data does not match with channel numbers.');
    end
    
    origEEG=EEG;
    
    % Remove extra events
    %evts_to_keep=find(strcmp('LBBS',{EEG.event.type}) | strcmp('LBSE',{EEG.event.type}) | strcmp('LBOB',{EEG.event.type}) | strcmp('LOBS',{EEG.event.type}) | strcmp('FTGO',{EEG.event.type}) | strcmp('LBEX',{EEG.event.type}) | strcmp('LEXT',{EEG.event.type}) | strcmp('FTGE',{EEG.event.type}));
    %EEG.event=EEG.event(evts_to_keep);
    %EEG=eeg_checkset(EEG,'eventconsistency');
    
    %% Delete discontinuous data from the raw data file
    % remove data after last task event
    latency=EEG.event(end).latency;
    % remove everything 1.5 seconds after the last event
    EEG = eeg_eegrej( EEG, [(latency+(1.5*EEG.srate)) EEG.pnts] );
    EEG = eeg_checkset( EEG );
    
    latency=EEG.event(1).latency;
    % remove everything until 1.5 seconds before the first event
    EEG = eeg_eegrej( EEG, [1 (latency-(1.5*EEG.srate))] );
    EEG = eeg_checkset( EEG );
        
    % Plot channel layout
    fig=figure();
    topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    saveas(fig, fullfile(subject_output_data_dir,'01-initial_ch_locations.png'));
        
    %% Filter data
    % Calculate filter order using the formula: m = dF / (df / fs), where m = filter order,
    % df = transition band width, dF = normalized transition width, fs = sampling rate
    % dF is specific for the window type. Hamming window dF = 3.3
    
    high_transband = highpass; % high pass transition band
    low_transband = 10; % low pass transition band
    
    hp_fl_order = 3.3 / (high_transband / EEG.srate);
    lp_fl_order = 3.3 / (low_transband / EEG.srate);
    
    % Round filter order to next higher even integer. Filter order is always even integer.
    if mod(floor(hp_fl_order),2) == 0
        hp_fl_order=floor(hp_fl_order);
    elseif mod(floor(hp_fl_order),2) == 1
        hp_fl_order=floor(hp_fl_order)+1;
    end
    
    if mod(floor(lp_fl_order),2) == 0
        lp_fl_order=floor(lp_fl_order)+2;
    elseif mod(floor(lp_fl_order),2) == 1
        lp_fl_order=floor(lp_fl_order)+1;
    end
    
    % Calculate cutoff frequency
    high_cutoff = highpass/2;
    low_cutoff = lowpass + (low_transband/2);
    
    % Performing high pass filtering
    EEG = eeg_checkset( EEG );
    EEG = pop_firws(EEG, 'fcutoff', high_cutoff, 'ftype', 'highpass',...
        'wtype', 'hamming', 'forder', hp_fl_order, 'minphase', 0);
    EEG = eeg_checkset( EEG );
    
    % Performing low pass filtering
    EEG = eeg_checkset( EEG );
    EEG = pop_firws(EEG, 'fcutoff', low_cutoff, 'ftype', 'lowpass',...
        'wtype', 'hamming', 'forder', lp_fl_order, 'minphase', 0);
    EEG = eeg_checkset( EEG );
    
    % Plot PSD
    fig=compute_and_plot_psd(EEG,1:EEG.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'05-filtered_psd.png'));
    
    %% NEAR Bad Channel Detection     
    try
        [EEG, flat_ch, lof_ch, periodo_ch, LOF_vec, thresh_lof_update] = NEAR_getBadChannels(EEG, isFlat, flatWin, isLOF, thresh_lof, dist_metric, isAdapt, ...
            isPeriodogram, frange, winsize, winov, pthresh, 0);
        save(fullfile(subject_output_data_dir, 'LOF_Values.mat'), 'LOF_vec'); % save .mat format
        disp('Bad Channel Detection is performed successfully');
        badChans = sort(unique(union(union(flat_ch, lof_ch),periodo_ch)));
    catch
        badChans=[];
    end

    if(~isempty(badChans))
        if(size(badChans,1) ~= 1)
            badChans = badChans';
        end
    end

    EEG = pop_select(EEG, 'nochannel', badChans);

    lof_flat_channels{s_idx}='';
    if numel(flat_ch)>0
        lof_flat_channels(s_idx)=join(cellfun(@(x) num2str(x(1)), num2cell(flat_ch,3), 'UniformOutput', false)',',');
    end
    lof_channels{s_idx}='';
    if numel(lof_ch)>0
        lof_channels(s_idx)=join(cellfun(@(x) num2str(x(1)), num2cell(lof_ch,3), 'UniformOutput', false)',',');
    end
    lof_periodo_channels{s_idx}='';
    if numel(periodo_ch)>0
        lof_periodo_channels(s_idx)=join(cellfun(@(x) num2str(x(1)), num2cell(periodo_ch,3), 'UniformOutput', false)',',');
    end
    lof_bad_channels{s_idx}='';
    if numel(badChans)>0
        lof_bad_channels(s_idx)=join(cellfun(@(x) num2str(x(1)), num2cell(badChans,3), 'UniformOutput', false)',',');
    end
    
    %% Save data after running filter and LOF function
    EEG = eeg_checkset( EEG );
    EEG = pop_editset(EEG, 'setname', strrep(data_file_name, '.vhdr', '_filtered_data'));
    EEG = pop_saveset( EEG,'filename',strrep(data_file_name, '.vhdr', '_filtered_data.set'),...
        'filepath', [subject_output_data_dir filesep '01_filtered_data' filesep]); % save .set format
    
    fig=figure();
    topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    saveas(fig, fullfile(subject_output_data_dir,'06-lof_removed.png'));
    
    %% Bad epochs correction/removal using ASR
    EEG_copy = EEG;
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off', ...
        'Highpass','off','BurstCriterion',rej_cutoff,'WindowCriterion',add_reject,'BurstRejection',rej_mode,'Distance','Euclidian');

    if(strcmp(rej_mode, 'on'))
        modified_mask = ~EEG.etc.clean_sample_mask;
    else
        modified_mask = sum(abs(EEG_copy.data-EEG.data),1) > 1e-10;
    end

    tot_samples_modified = (length(find(modified_mask)) * 100) / EEG_copy.pnts;
    tot_samples_modified = round(tot_samples_modified * 100) / 100;
    asr_tot_samples_modified(s_idx)=tot_samples_modified;
    change_in_RMS = -(mean(rms(EEG.data,2)) - mean(rms(EEG_copy.data,2))*100)/mean(rms(EEG_copy.data,2)); % in percentage
    change_in_RMS = round(change_in_RMS * 100) / 100;
    asr_change_in_RMS(s_idx) =change_in_RMS;
    fprintf('\nArtifacted epochs are corrected by ASR algorithm\n');
    
    %% Save data after running ASR function
    EEG = eeg_checkset( EEG );
    EEG = pop_editset(EEG, 'setname', strrep(data_file_name, '.vhdr', '_asr_data'));
    EEG = pop_saveset( EEG,'filename',strrep(data_file_name, '.vhdr', '_asr_data.set'),...
        'filepath', [subject_output_data_dir filesep '02_near_data' filesep]); % save .set format

    fig=compute_and_plot_psd(EEG,1:EEG.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'07-asr_psd.png'));
    
    %% STEP 8: Prepare data for ICA
    EEG_copy=EEG; % make a copy of the dataset
    EEG_copy = eeg_checkset(EEG_copy);
    
    % Perform 1Hz high pass filter on copied dataset
    transband = 1;
    fl_cutoff = transband/2;
    fl_order = 3.3 / (transband / EEG.srate);
    
    if mod(floor(fl_order),2) == 0
        fl_order=floor(fl_order);
    elseif mod(floor(fl_order),2) == 1
        fl_order=floor(fl_order)+1;
    end
    
    EEG_copy = pop_firws(EEG_copy, 'fcutoff', fl_cutoff,...
        'ftype', 'highpass', 'wtype', 'hamming', 'forder', fl_order,...
        'minphase', 0);
    EEG_copy = eeg_checkset(EEG_copy);
    
    % Create 1 second epoch
    % insert temporary marker 1 second apart and create epochs
    EEG_copy=eeg_regepochs(EEG_copy,'recurrence', 1, 'limits',[0 1],...
        'rmbase', [NaN], 'eventtype', '999'); 
    EEG_copy = eeg_checkset(EEG_copy);
    
    % Find bad epochs and delete them from dataset
    % [lower upper] threshold limit(s) in mV.
    vol_thrs = [-1000 1000]; 
    
    % Find channel/s with xx% of artifacted 1-second epochs and delete them
    chanCounter = 1; ica_prep_badChans = [];
    numEpochs =EEG_copy.trials; % find the number of epochs
    all_bad_channels=0;
    
    for ch=1:EEG_copy.nbchan
        % Find artifaceted epochs by detecting outlier voltage
        EEG_copy = pop_eegthresh(EEG_copy,1, ch, vol_thrs(1), vol_thrs(2),...
            EEG_copy.xmin, EEG_copy.xmax, 0, 0);
        EEG_copy = eeg_checkset( EEG_copy );
        
        % Find number of artifacted epochs
        EEG_copy = eeg_checkset( EEG_copy );
        EEG_copy = eeg_rejsuperpose( EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
        artifacted_epochs=EEG_copy.reject.rejglobal;
        
        % Find bad channel / channel with more than 20% artifacted epochs
        if sum(artifacted_epochs) > (numEpochs*20/100)
            ica_prep_badChans(chanCounter) = ch;
            chanCounter=chanCounter+1;
        end
    end
    
    % If all channels are bad, save the dataset at this stage and ignore the remaining of the preprocessing.
    if numel(ica_prep_badChans)==EEG.nbchan || numel(ica_prep_badChans)+1==EEG.nbchan
        all_bad_channels=1;
        warning(['No usable data for datafile', data_file_name]);        
    else
        % Reject bad channel - channel with more than xx% artifacted epochs
        EEG_copy = pop_select( EEG_copy,'nochannel', ica_prep_badChans);
        EEG_copy = eeg_checkset(EEG_copy);
    end
    
    if numel(ica_prep_badChans)==0
        ica_preparation_bad_channels{s_idx}='0';
    else
        ica_preparation_bad_channels{s_idx}=num2str(ica_prep_badChans);
    end
    
    if all_bad_channels == 1
        length_ica_data(s_idx)=0;
        total_ICs(s_idx)=0;
        ICs_removed{s_idx}='0';
        total_epochs_before_artifact_rejection(s_idx)=0;
        total_epochs_after_artifact_rejection(s_idx)=0;
        total_channels_interpolated(s_idx)=0;
        continue % ignore rest of the processing and go to next datafile
    end
    
    % Find the artifacted epochs across all channels and reject them before doing ICA.
    EEG_copy = pop_eegthresh(EEG_copy,1, 1:EEG_copy.nbchan, vol_thrs(1),...
        vol_thrs(2), EEG_copy.xmin, EEG_copy.xmax,0,0);
    EEG_copy = eeg_checkset(EEG_copy);
    
    % Find the number of artifacted epochs and reject them
    EEG_copy = eeg_checkset(EEG_copy);
    EEG_copy = eeg_rejsuperpose(EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
    reject_artifacted_epochs=EEG_copy.reject.rejglobal;
    EEG_copy = pop_rejepoch(EEG_copy, reject_artifacted_epochs, 0);
    
    fig=compute_and_plot_psd(EEG_copy, 1:EEG_copy.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'08-ica_copy_epochs_psd.png'));    
    
    %% Run ICA
    length_ica_data(s_idx)=EEG_copy.trials; % length of data (in second) fed into ICA
    EEG_copy = eeg_checkset(EEG_copy);
    EEG_copy = pop_runica(EEG_copy, 'icatype', 'runica', 'extended', 1,...
        'stop', 1E-7, 'interupt','off');
%     EEG_copy = pop_runica(EEG_copy, 'icatype', 'sobi');
    
    EEG_copy = eeg_checkset(EEG_copy);
    EEG_copy = pop_editset(EEG_copy, 'setname',  strrep(data_file_name, '.vhdr', '_ica'));
    EEG_copy = pop_saveset(EEG_copy, 'filename', strrep(data_file_name, '.vhdr', '_ica.set'),...
        'filepath', [subject_output_data_dir filesep '03_ica_data' filesep ]); % save .set format
    
    % Find the ICA weights that would be transferred to the original dataset
    ICA_WINV=EEG_copy.icawinv;
    ICA_SPHERE=EEG_copy.icasphere;
    ICA_WEIGHTS=EEG_copy.icaweights;
    ICA_CHANSIND=EEG_copy.icachansind;
    
    % If channels were removed from copied dataset during preparation of ica, then remove
    % those channels from original dataset as well before transferring ica weights.
    EEG = eeg_checkset(EEG);
    EEG = pop_select(EEG,'nochannel', ica_prep_badChans);
    
    % Transfer the ICA weights of the copied dataset to the original dataset
    EEG.icawinv=ICA_WINV;
    EEG.icasphere=ICA_SPHERE;
    EEG.icaweights=ICA_WEIGHTS;
    EEG.icachansind=ICA_CHANSIND;
    EEG = eeg_checkset(EEG);
    
    %% Run adjusted-adjust to find artifacted ICA components
    badICs=[];    
    if size(EEG_copy.icaweights,1) == size(EEG_copy.icaweights,2)
        figure()
        badICs = adjusted_ADJUST(EEG_copy, [[subject_output_data_dir filesep '03_ica_data' filesep] strrep(data_file_name, '.vhdr', '_adjust_report')]);        
        close all;
    else % if rank is less than the number of electrodes, throw a warning message
        warning('The rank is less than the number of electrodes. ADJUST will be skipped. Artefacted ICs will have to be manually rejected for this participant');
    end
    adjadj_fig_fname=sprintf('%s.jpg', subject);
    if exist(adjadj_fig_fname,'file')==2
        movefile(adjadj_fig_fname, fullfile(subject_output_data_dir, 'adjusted_adjust.jpg'));
    end
    
    % Mark the bad ICs found by ADJUST
    for ic=1:length(badICs)
        EEG.reject.gcompreject(1, badICs(ic))=1;
        EEG = eeg_checkset(EEG);
    end
    total_ICs(s_idx)=size(EEG.icasphere, 1);
    if numel(badICs)==0
        ICs_removed{s_idx}='0';
    else
        ICs_removed{s_idx}=num2str(double(badICs));
    end    
    
    %% Save dataset after ICA
    EEG = eeg_checkset(EEG);
    EEG = pop_editset(EEG, 'setname',  strrep(data_file_name, '.vhdr', '_ica_data'));
    EEG = pop_saveset(EEG, 'filename', strrep(data_file_name, '.vhdr', '_ica_data.set'),...
        'filepath', [subject_output_data_dir filesep '03_ica_data' filesep ]); % save .set format
    
    %% Remove artifacted ICA components from data
    all_bad_ICs=0;
    ICs2remove=find(EEG.reject.gcompreject); % find ICs to remove
    
    % If all ICs and bad, save data at this stage and ignore rest of the preprocessing for this subject.
    if numel(ICs2remove)==total_ICs(s_idx)
        all_bad_ICs=1;
        warning(['No usable data for datafile', data_file_name]);        
    else
        EEG = eeg_checkset( EEG );
        EEG = pop_subcomp( EEG, ICs2remove, 0); % remove ICs from dataset
    end
    
    if all_bad_ICs==1
        total_epochs_before_artifact_rejection(s_idx)=0;
        total_epochs_after_artifact_rejection(s_idx)=0;
        total_channels_interpolated(s_idx)=0;
        continue % ignore rest of the processing and go to next datafile
    end
    
    fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'09-ica_art_rej_psd.png'));
    
    
    
    %% Segment data into fixed length epochs
    EEG = eeg_checkset(EEG);
    EEG1 = pop_epoch(EEG, event_markers_1, epoch_length_1, 'epochinfo', 'yes');
    
    % Delete non-time-locking events from epoched data (08/16/2020 updated)
    for epochIdx = 1:length(EEG1.epoch)
        allEventIdx    = 1:length(EEG1.epoch(epochIdx).event);
        if length(allEventIdx) == 1 % If there is only 1 event, it must have latency zero.
            continue
        else
            zeroLatencyIdx = find(cell2mat(EEG1.epoch(epochIdx).eventlatency) == 0);
            EEG1.epoch(epochIdx).event         = EEG1.epoch(epochIdx).event(zeroLatencyIdx);
            EEG1.epoch(epochIdx).eventlatency  = {0};
            EEG1.epoch(epochIdx).eventduration = 0;
            EEG1.epoch(epochIdx).eventchannel = EEG1.epoch(epochIdx).eventchannel(zeroLatencyIdx);
            EEG1.epoch(epochIdx).eventbvtime = EEG1.epoch(epochIdx).eventbvtime(zeroLatencyIdx);
            EEG1.epoch(epochIdx).eventbvmknum = EEG1.epoch(epochIdx).eventbvmknum(zeroLatencyIdx);
            EEG1.epoch(epochIdx).eventvisible = EEG1.epoch(epochIdx).eventvisible(zeroLatencyIdx);
            EEG1.epoch(epochIdx).eventtype     = EEG1.epoch(epochIdx).eventtype(zeroLatencyIdx);
            EEG1.epoch(epochIdx).eventcode     = EEG1.epoch(epochIdx).eventcode(zeroLatencyIdx);
            EEG1.epoch(epochIdx).eventurevent  = EEG1.epoch(epochIdx).eventurevent(zeroLatencyIdx);
            EEG1.epoch(epochIdx).eventattention = EEG1.epoch(epochIdx).eventattention(zeroLatencyIdx);
        end
    end
    validEventIdx  = [EEG1.epoch.event];
    deleteEventIdx = setdiff(1:length(EEG1.event), validEventIdx);
    EEG1.event(deleteEventIdx) = [];
    for epochIdx = 1:length(EEG1.epoch)
        EEG1.epoch(epochIdx).event = epochIdx;
    end

    fig=compute_and_plot_psd(EEG1, 1:EEG1.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'10-epoch_psd_eye.png'));
    
    total_eye_epochs_before_artifact_rejection(s_idx)=EEG1.trials;
    
    %% Artifact rejection
    all_bad_epochs=0;
    chans=[]; chansidx=[];chans_labels2=[];
    chans_labels2=cell(1,EEG1.nbchan);
    for i=1:EEG1.nbchan
        chans_labels2{i}= EEG1.chanlocs(i).labels;
    end
    [chans,chansidx] = ismember(frontal_channels, chans_labels2);
    frontal_channels_idx = chansidx(chansidx ~= 0);
    badChans = zeros(EEG1.nbchan, EEG1.trials);
    badepoch=zeros(1, EEG1.trials);
    
    % find artifaceted epochs by detecting outlier voltage in the specified channels list and remove epoch if artifacted in those channels
    for ch =1:length(frontal_channels_idx)
        EEG1 = pop_eegthresh(EEG1,1, frontal_channels_idx(ch), volt_threshold(1), volt_threshold(2), EEG1.xmin, EEG1.xmax,0,0);
        EEG1 = eeg_checkset( EEG1 );
        EEG1 = eeg_rejsuperpose( EEG1, 1, 1, 1, 1, 1, 1, 1, 1);
        badChans(ch,:) = EEG1.reject.rejglobal;
    end
    for ii=1:size(badChans, 2)
        badepoch(ii)=sum(badChans(:,ii));
    end
    badepoch=logical(badepoch);
    
    % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
    if sum(badepoch)==EEG1.trials || sum(badepoch)+1==EEG1.trials
        all_bad_epochs=1;
        warning(['No usable eye data for datafile', data_file_name]);                
    else
        EEG1 = pop_rejepoch( EEG1, badepoch, 0);
        EEG1 = eeg_checkset(EEG1);
    end

    if all_bad_epochs==0
        % Interpolate artifacted data for all reaming channels
        badChans = zeros(EEG1.nbchan, EEG1.trials);
        % Find artifacted epochs by detecting outlier voltage but don't remove
        for ch=1:EEG1.nbchan
            EEG1 = pop_eegthresh(EEG1,1, ch, volt_threshold(1), volt_threshold(2), EEG1.xmin, EEG1.xmax,0,0);
            EEG1 = eeg_checkset(EEG1);
            EEG1 = eeg_rejsuperpose(EEG1, 1, 1, 1, 1, 1, 1, 1, 1);
            badChans(ch,:) = EEG1.reject.rejglobal;
        end
        tmpData = zeros(EEG1.nbchan, EEG1.pnts, EEG1.trials);
        for e = 1:EEG1.trials
            % Initialize variables EEGe and EEGe_interp;
            EEGe1 = []; EEGe1_interp = []; badChanNum = [];
            % Select only this epoch (e)
            EEGe1 = pop_selectevent( EEG1, 'epoch', e, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
            badChanNum = find(badChans(:,e)==1); % find which channels are bad for this epoch
            if length(badChanNum) < round((10/100)*EEG1.nbchan)% check if more than 10% are bad
                EEGe1_interp = eeg_interp(EEGe1,badChanNum); %interpolate the bad channels for this epoch
                tmpData(:,:,e) = EEGe1_interp.data; % store interpolated data into matrix
            end
        end
        EEG1.data = tmpData; % now that all of the epochs have been interpolated, write the data back to the main file

        % If more than 10% of channels in an epoch were interpolated, reject that epoch
        badepoch=zeros(1, EEG1.trials);
        for ei=1:EEG1.trials
            NumbadChan = length(find(badChans(:,ei)==1)); % find how many channels are bad in an epoch
            if sum(NumbadChan) >= round((10/100)*EEG1.nbchan)% check if more than 10% are bad
                badepoch (ei)= sum(NumbadChan);
            end
        end
        badepoch=logical(badepoch);
    end
    
    % Find indices of events where attention is 'not_attentive'
    idx_remove = find(strcmp({EEG1.event.attention}, 'not_attentive'));
    badepoch(idx_remove)=1;

    % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
    if sum(badepoch)==EEG1.trials || sum(badepoch)+1==EEG1.trials
        all_bad_epochs=1;
        warning(['No usable eye data for datafile', data_file_name]);                
    else
        EEG1 = pop_rejepoch(EEG1, badepoch, 0);
        EEG1 = eeg_checkset(EEG1);
    end
    
    % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
    if sum(EEG1.reject.rejthresh)==EEG1.trials || sum(EEG1.reject.rejthresh)+1==EEG1.trials
        all_bad_epochs=1;
        warning(['No usable eye data for datafile', data_file_name]);            
    else
        EEG1 = pop_rejepoch(EEG1,(EEG1.reject.rejthresh), 0);
        EEG1 = eeg_checkset(EEG1);
    end
    
    % if all epochs are found bad during artifact rejection
    if all_bad_epochs==1
        total_eye_epochs_after_artifact_rejection(s_idx)=0;
        total_eye_channels_interpolated(s_idx)=0;
    else
        total_eye_epochs_after_artifact_rejection(s_idx)=EEG1.trials;        
    
        %% Interpolation  
        total_eye_channels_interpolated(s_idx)=length(origEEG.chanlocs)-length(EEG1.chanlocs);
        EEG1 = pop_interp(EEG1, origEEG.chanlocs, interp_type);
        fprintf('\nMissed channels are spherically interpolated\n');

        %% Re-referencing
        EEG1 = pop_reref( EEG1, []);

        fig=compute_and_plot_psd(EEG1, 1:EEG1.nbchan);
        saveas(fig, fullfile(subject_output_data_dir,'11-art_rej_reref_psd_eye.png'));    

        %% Save processed data
        EEG1 = eeg_checkset(EEG1);
        EEG1 = pop_editset(EEG1, 'setname',  strrep(data_file_name, '.vhdr', '_rereferenced_data_eye'));
        EEG1 = pop_saveset(EEG1, 'filename', strrep(data_file_name, '.vhdr', '_rereferenced_data_eye.set'),...
            'filepath', [subject_output_data_dir filesep '04_rereferenced_data' filesep ]); % save .set format
    end
    
    
    
    
    EEG2 = pop_epoch(EEG, event_markers_2, epoch_length_2, 'epochinfo', 'yes');
    
    % Delete non-time-locking events from epoched data (08/16/2020 updated)
    for epochIdx = 1:length(EEG2.epoch)
        allEventIdx    = 1:length(EEG2.epoch(epochIdx).event);
        if length(allEventIdx) == 1 % If there is only 1 event, it must have latency zero.
            continue
        else
            zeroLatencyIdx = find(cell2mat(EEG2.epoch(epochIdx).eventlatency) == 0);
            EEG2.epoch(epochIdx).event         = EEG2.epoch(epochIdx).event(zeroLatencyIdx);
            EEG2.epoch(epochIdx).eventtype     = EEG2.epoch(epochIdx).eventtype(zeroLatencyIdx);
            EEG2.epoch(epochIdx).eventlatency  = {0};
            EEG2.epoch(epochIdx).eventattention = EEG2.epoch(epochIdx).eventattention(zeroLatencyIdx);
            EEG2.epoch(epochIdx).eventduration = 0;
            EEG2.epoch(epochIdx).eventchannel = EEG2.epoch(epochIdx).eventchannel(zeroLatencyIdx);
            EEG2.epoch(epochIdx).eventbvtime = EEG2.epoch(epochIdx).eventbvtime(zeroLatencyIdx);
            EEG2.epoch(epochIdx).eventbvmknum = EEG2.epoch(epochIdx).eventbvmknum(zeroLatencyIdx);
            EEG2.epoch(epochIdx).eventvisible = EEG2.epoch(epochIdx).eventvisible(zeroLatencyIdx);
            EEG2.epoch(epochIdx).eventcode     = EEG2.epoch(epochIdx).eventcode(zeroLatencyIdx);
            EEG2.epoch(epochIdx).eventurevent  = EEG2.epoch(epochIdx).eventurevent(zeroLatencyIdx);
        end
    end
    validEventIdx  = [EEG2.epoch.event];
    deleteEventIdx = setdiff(1:length(EEG2.event), validEventIdx);
    EEG2.event(deleteEventIdx) = [];
    for epochIdx = 1:length(EEG2.epoch)
        EEG2.epoch(epochIdx).event = epochIdx;
    end

    fig=compute_and_plot_psd(EEG2, 1:EEG2.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'10-epoch_psd_task.png'));
    
    total_task_epochs_before_artifact_rejection(s_idx)=EEG2.trials;
    
    %% Artifact rejection
    all_bad_epochs=0;
    chans=[]; chansidx=[];chans_labels2=[];
    chans_labels2=cell(1,EEG2.nbchan);
    for i=1:EEG2.nbchan
        chans_labels2{i}= EEG2.chanlocs(i).labels;
    end
    [chans,chansidx] = ismember(frontal_channels, chans_labels2);
    frontal_channels_idx = chansidx(chansidx ~= 0);
    badChans = zeros(EEG2.nbchan, EEG2.trials);
    badepoch=zeros(1, EEG2.trials);
    
    % find artifaceted epochs by detecting outlier voltage in the specified channels list and remove epoch if artifacted in those channels
    for ch =1:length(frontal_channels_idx)
        EEG2 = pop_eegthresh(EEG2,1, frontal_channels_idx(ch), volt_threshold(1), volt_threshold(2), EEG2.xmin, EEG2.xmax,0,0);
        EEG2 = eeg_checkset( EEG2 );
        EEG2 = eeg_rejsuperpose( EEG2, 1, 1, 1, 1, 1, 1, 1, 1);
        badChans(ch,:) = EEG2.reject.rejglobal;
    end
    for ii=1:size(badChans, 2)
        badepoch(ii)=sum(badChans(:,ii));
    end
    badepoch=logical(badepoch);
        
    
    % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
    if sum(badepoch)==EEG2.trials || sum(badepoch)+1==EEG2.trials
        all_bad_epochs=1;
        warning(['No usable task data for datafile', data_file_name]);                
    else
        EEG2 = pop_rejepoch( EEG2, badepoch, 0);
        EEG2 = eeg_checkset(EEG2);
    end

    if all_bad_epochs==0
        % Interpolate artifacted data for all reaming channels
        badChans = zeros(EEG2.nbchan, EEG2.trials);
        % Find artifacted epochs by detecting outlier voltage but don't remove
        for ch=1:EEG2.nbchan
            EEG2 = pop_eegthresh(EEG2,1, ch, volt_threshold(1), volt_threshold(2), EEG2.xmin, EEG2.xmax,0,0);
            EEG2 = eeg_checkset(EEG2);
            EEG2 = eeg_rejsuperpose(EEG2, 1, 1, 1, 1, 1, 1, 1, 1);
            badChans(ch,:) = EEG2.reject.rejglobal;
        end
        tmpData = zeros(EEG2.nbchan, EEG2.pnts, EEG2.trials);
        for e = 1:EEG2.trials
            % Initialize variables EEGe and EEGe_interp;
            EEGe2 = []; EEGe2_interp = []; badChanNum = [];
            % Select only this epoch (e)
            EEGe2 = pop_selectevent( EEG2, 'epoch', e, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
            badChanNum = find(badChans(:,e)==1); % find which channels are bad for this epoch
            if length(badChanNum) < round((10/100)*EEG2.nbchan)% check if more than 10% are bad
                EEGe2_interp = eeg_interp(EEGe2,badChanNum); %interpolate the bad channels for this epoch
                tmpData(:,:,e) = EEGe2_interp.data; % store interpolated data into matrix
            end
        end
        EEG2.data = tmpData; % now that all of the epochs have been interpolated, write the data back to the main file

        % If more than 10% of channels in an epoch were interpolated, reject that epoch
        badepoch=zeros(1, EEG2.trials);
        for ei=1:EEG2.trials
            NumbadChan = length(find(badChans(:,ei)==1)); % find how many channels are bad in an epoch
            if sum(NumbadChan) >= round((10/100)*EEG2.nbchan)% check if more than 10% are bad
                badepoch (ei)= sum(NumbadChan);
            end
        end
        badepoch=logical(badepoch);
    end
    
    % Find indices of events where attention is 'not_attentive'
    idx_remove = find(strcmp({EEG2.event.attention}, 'not_attentive'));
    badepoch(idx_remove)=1;

    % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
    if sum(badepoch)==EEG2.trials || sum(badepoch)+1==EEG2.trials
        all_bad_epochs=1;
        warning(['No usable task data for datafile', data_file_name]);                
    else
        EEG2 = pop_rejepoch(EEG2, badepoch, 0);
        EEG2 = eeg_checkset(EEG2);
    end
    
    % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
    if sum(EEG2.reject.rejthresh)==EEG2.trials || sum(EEG2.reject.rejthresh)+1==EEG2.trials
        all_bad_epochs=1;
        warning(['No usable task data for datafile', data_file_name]);            
    else
        EEG2 = pop_rejepoch(EEG2,(EEG2.reject.rejthresh), 0);
        EEG2 = eeg_checkset(EEG2);
    end
    
    % if all epochs are found bad during artifact rejection
    if all_bad_epochs==1
        total_task_epochs_after_artifact_rejection(s_idx)=0;
        total_task_channels_interpolated(s_idx)=0;
    else
        total_task_epochs_after_artifact_rejection(s_idx)=EEG2.trials;        
    
        %% Interpolation  
        total_task_channels_interpolated(s_idx)=length(origEEG.chanlocs)-length(EEG2.chanlocs);
        EEG2 = pop_interp(EEG2, origEEG.chanlocs, interp_type);
        fprintf('\nMissed channels are spherically interpolated\n');

        %% Re-referencing
        EEG2 = pop_reref( EEG2, []);

        fig=compute_and_plot_psd(EEG2, 1:EEG2.nbchan);
        saveas(fig, fullfile(subject_output_data_dir,'11-art_rej_reref_psd_task.png'));    

        %% Save processed data
        EEG2 = eeg_checkset(EEG2);
        EEG2 = pop_editset(EEG2, 'setname',  strrep(data_file_name, '.vhdr', '_rereferenced_data_task'));
        EEG2 = pop_saveset(EEG2, 'filename', strrep(data_file_name, '.vhdr', '_rereferenced_data_task.set'),...
            'filepath', [subject_output_data_dir filesep '04_rereferenced_data' filesep ]); % save .set format
    end
    
    close all;
end

%% Create the report table for all the data files with relevant preprocessing outputs.
report_table=table(study_info.participant_info.participant_id,...
    lof_flat_channels', lof_channels', lof_periodo_channels', lof_bad_channels',...
    asr_tot_samples_modified', asr_change_in_RMS', ica_preparation_bad_channels',...
    length_ica_data', total_ICs', ICs_removed', total_eye_epochs_before_artifact_rejection',...
    total_eye_epochs_after_artifact_rejection',total_eye_channels_interpolated',...
    total_task_epochs_before_artifact_rejection', total_task_epochs_after_artifact_rejection',...
    total_task_channels_interpolated');

report_table.Properties.VariableNames={'subject','lof_flat_channels', 'lof_channels', ...
    'lof_periodo_channels', 'lof_bad_channels', 'asr_tot_samples_modified', 'asr_change_in_RMS',...
    'ica_preparation_bad_channels', 'length_ica_data', 'total_ICs', 'ICs_removed',...
    'total_eye_epochs_before_artifact_rejection', 'total_eye_epochs_after_artifact_rejection',...
    'total_eye_channels_interpolated', 'total_task_epochs_before_artifact_rejection',...
    'total_task_epochs_after_artifact_rejection', 'total_task_channels_interpolated'};
writetable(report_table, fullfile(study_info.data_dir, 'derivatives', 'NEARICA_behav_v3', ['NEARICA_preprocessing_report_', datestr(now,'dd-mm-yyyy'),'.csv']));
%end