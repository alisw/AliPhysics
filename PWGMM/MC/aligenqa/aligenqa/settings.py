# A list of the estimators which will be included in the plots (if available in the input file)
considered_ests = [
    'EtaLt05',
    'EtaLt08',
    'EtaLt15',
    'Eta08_15',
    'V0M',
    # 'V0A',
    'V0C',
    'ZDC',
    'nMPI',
    'Q2',
    'spherocity',
    'sphericity'
]

# Considered triggers
considered_triggers = [
    'Inel',
    'InelGt0',
    'V0AND'
]

# Ranges of percentiles which should be considered in the plots These
# ranges are then translated into ranges of bins in multiplicity, nMPI
# or whatever is applicable
std_perc_bins = [(1, 0.7), (.5, .4), (.1, .05), (0.001, 0.0)]
percentile_bins = {
    'EtaLt05': std_perc_bins,
    'EtaLt08': std_perc_bins,
    'EtaLt15': std_perc_bins,
    'Eta08_15': std_perc_bins,
    'V0M': std_perc_bins,
    'V0A': std_perc_bins,
    'V0C': std_perc_bins,
    'ZDC': [(1, 0.7), (.7, .3), (.3, .05), (0.001, 0.0)],
    'nMPI': [(1, 0.7), (.7, .4), (.3, .05), (0.001, 0.0)],
    'Q2': [(1, 0.7), (.7, .4), (.3, .05), (0.001, 0.0)],
    'spherocity': [(1, 0.7), (.7, .4), (.3, .05), (0.001, 0.0)],
    'sphericity': [(1, 0.7), (.7, .4), (.3, .05), (0.001, 0.0)],
}
