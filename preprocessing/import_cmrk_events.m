function EEG = import_cmrk_events(EEG, cmrk_filepath)

fid = fopen(cmrk_filepath, 'r');
if fid == -1
    error('Cannot open .cmrk file: %s', cmrk_filepath);
end

% Build an index of EEG events based on type and latency for quick matching
EEG_event_idx = containers.Map();
for idx = 1:numel(EEG.event)
    key = sprintf('%s_%d', EEG.event(idx).type, EEG.event(idx).latency);
    EEG_event_idx(key) = idx;
    EEG.event(idx).attention = 'not_applicable'; % default initialization
end

% Process cmrk file
while ~feof(fid)
    line = fgetl(fid);

    if startsWith(line, 'Mk')
        tokens = regexp(line, 'Mk\d+=Event,([^,]+),(\d+),\d+,?(\d*)', 'tokens');
        if ~isempty(tokens)
            tokens = tokens{1};
            event_type = strtrim(tokens{1});
            event_latency = str2double(tokens{2});
            key = sprintf('%s_%d', event_type, event_latency);

            if isKey(EEG_event_idx, key)
                EEG_idx = EEG_event_idx(key);
                EEG.event(EEG_idx).attention = 'attentive';
                if length(tokens) >= 3 && ~isempty(tokens{3})
                    att_val = str2double(tokens{3});
                    if att_val == 0
                        EEG.event(EEG_idx).attention = 'not_attentive';                    
                    end
                end
            else
                warning('Event type %s at latency %d in cmrk not found in EEG events.', event_type, event_latency);
            end
        end
    end
end

fclose(fid);

EEG = eeg_checkset(EEG, 'eventconsistency');

end
