import sys
import os.path as op
import mne
import pandas as pd
import scipy
import numpy as np

from zapline_iter import zapline_until_gone

try:
    base_dir = sys.argv[1]
    pipeline = sys.argv[2]
    step_path = sys.argv[3]
except:
    print("incorrect arguments")
    sys.exit()

subjects=pd.read_csv(op.join(base_dir, 'data', 'participants.tsv'), sep='\t')
subject_ids = subjects['participant_id']

for subject_id in subject_ids:

    eeg_path = op.join(base_dir,'data/derivatives',pipeline,subject_id, step_path)
    out_path = op.join(base_dir, 'data/derivatives', pipeline, subject_id, '05_zapped_data')

    f_name=op.join(eeg_path,'%s_rereferenced_data_eye.set' % subject_id)
    if op.exists(f_name):
        epochs=mne.read_epochs_eeglab(f_name)
        [data, iterations]=zapline_until_gone(epochs.get_data(), 50, epochs.info['sfreq'], viz=True,
                                              prefix=op.join(base_dir, 'data/derivatives', pipeline, subject_id,
                                                             "zapline_iter_eye"),
                                              max_iter=3)

        scipy.io.savemat(op.join(out_path,'zapline_eye.mat'),{'data':data,'iteration':iterations})

    f_name = op.join(eeg_path, '%s_rereferenced_data_task.set' % subject_id)
    if op.exists(f_name):
        epochs = mne.read_epochs_eeglab(f_name)
        [data, iterations] = zapline_until_gone(epochs.get_data(), 50, epochs.info['sfreq'], viz=True,
                                                prefix=op.join(base_dir, 'data/derivatives', pipeline, subject_id,
                                                               "zapline_iter_task"),
                                                max_iter=3)

        scipy.io.savemat(op.join(out_path, 'zapline_task.mat'), {'data': data, 'iteration': iterations})
