# Mouse_CNT
Optogenetics and electrophysiology recordings of neural and muscle activities in the mouse.
Data and MATLAB scripts used to generate results presented in the JP paper are provided here.

### Spike-triggered averages of EMGs
Onset latency (`OnsetLatency_all.mat`) and peak width at half maximum (`PWHM_all.mat`) were used two-dimensional feature space of all spike triggered averages to characterize synchrony effects in our dataset.

### Preferred and Non-preferred pairs
`Similarity_MF_MS.xlsx` contains scalar product values of preferred pairs and non-preferred pairs between muscle fields and muscle synergies across all mice.

### Clustering of muscle synergies
`cluster_ms.mat` contains clusters of muscle synergies. There are multiple fields within the data structure: 
- `subject`: mouse id
- `property`: preferred (`p`) or unmatched (`u`)
- `mf_subjects`: mouse id of the muscle fields
- `mf`: muscle fields
- `data`: muscle synergies 
Clusters and plots can be generated from script `MatchMFintoWClusters.m`

### Pairwise similarity between vectors
The pairwise similarity results can be found in the data structure `SST_vs_SSAT.mat` for conditions in\between SST and SSAT and `SST_vs_CoST.mat` for conditions in\between SST and CoST. Specifically, within `SST_vs_SSAT.mat`, you will find field `SP_SST` and `SP_SSAT` representing the scalar product value of best-match pairs between muscle field and muscle synergies in SST and SSAT, respectively. Field `SP_Un_SST` and `SP_Un_SSAT` represent the scalar product value of unmatch pairs in SST and SSAT. Field `SP_SST_SSAT` contains scalar product values of muscle fields across SST and SSAT conditions. Field `SP_Syn_SST_SSAT` contains scalar product values of muscle synergies across SST and SSAT conditions.

### Temporal correlation between neural firing rate and synergy activities.
The inidividual mean firing rate and synergy activities can be found in `fr.mat` and `Coef.mat`, which could be used to generate the correlation results between individual firing rate and synergy activities of the same mouse. The synergies can be obtained from `mouse_19_12_syn_SS_T.mat`. The firing rates and synergy activities across all mice can be found in `FiringRate_all.mat` and `Ccorres.mat`, as well as `Cprefer.mat`. The script used to generate the above results is `CorrCampFR.m`
