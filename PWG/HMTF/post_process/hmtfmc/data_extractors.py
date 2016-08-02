"""
This file provides functions to extract data from the "Sums" histogram. The functions should, return (lists of)
plottables.
"""

from rootpy import asrootpy

from utils import \
    gen_random_name,\
    get_est_dirs,\
    make_estimator_title,\
    remap_x_values,\
    remove_zero_value_points,\
    remove_points_with_equal_x,\
    remove_points_with_x_err_gt_1NchRef,\
    remove_non_mutual_points

import ROOT


def get_meanpt_vs_estmult(resutlts_est_dir, pids):
    """
    Create a 1Dprofile for the given pids and the given estimator name
    """
    # find the mult vs pt histograms for the given pids
    mult_vs_pts = []
    for pid in pids:
        mult_vs_pts.append(asrootpy(getattr(resutlts_est_dir.mult_pt, pid)))
    profx = sum(mult_vs_pts).ProfileX()
    profx.name = gen_random_name()
    return profx


def get_dNdeta_in_classifier_bin_interval(sums_classifier_dir, event_counter, classifier_bin_interval):
    """
    Get dN/deta for a given interval of classifier bin indices
    Parameters
    ----------
    sums_classifier_dir : TList
        Sums directory of a classifier
    event_counter : Hist1D
        Event counter histogram with the classifier value on the xaxis
    classifier_bin_interval : list
        classifier value bin edges given as bin indices
    Returns
    -------
    Hist1D
    """
    hist_name = "eta_classifier_{0}".format(sums_classifier_dir.GetName())
    h2d = asrootpy(sums_classifier_dir.FindObject(hist_name))
    if not h2d:
        raise ValueError("Could not find histogram {0}".format(hist_name))
    h2d.yaxis.set_range(classifier_bin_interval[0], classifier_bin_interval[1])
    h = asrootpy(h2d.projection_x(gen_random_name()))
    h.title = "{0} - {1} %".format(100 * classifier_bin_interval[0], 100 * classifier_bin_interval[1])
    # scale by the number of events in this mult_interval and bin width
    try:
        h.Scale((1.0 /
                 float(event_counter.Integral(classifier_bin_interval[0], classifier_bin_interval[1]))),
                "width")
    except ZeroDivisionError:
        # If this happens, we have empty bins in dN/deta! The stats must suck!
        raise ZeroDivisionError("Your statistics are terrible! Consider increasing the classifier value interval to avoid this")
    return h


def get_identified_vs_mult(h3d, pdg):
    """
    Return 1D counter histogram of identified particles vs N_ch^est
    Parameters
    ----------
    h3d: Hist3D
         x: est_mult; y: pT; z: pids1
    pdg: str
         pdg code as string
    Return
    ------
    Hist1D:
         x: Nch_est y: counts
    """
    pid_bin = h3d.zaxis.find_bin(pdg)
    if pid_bin == 0:
        raise ValueError("given pdg ({0}) does not exist in histogram".format(pdg))

    h3d.zaxis.SetRange(pid_bin, pid_bin)
    h = asrootpy(h3d.Project3D("yx").ProjectionX())
    h.SetName(gen_random_name())
    return h


def get_correlation_histogram(sums, classifier1, classifier2):
    """
    Returns the correlation histogram between classifier 1 and 2

    Parameters
    ----------
    sums : TList
           Sums directory
    classifier1 : str
           Name of estimator 1
    classifier2 : str
           Name of estimator 2
    Returns
    -------
    Hist2D :
        Hist2D with classifier1 on the x- and classifier2 on the y-axis_title

    Raises
    ------
    ValueError :
        If either of the classifiers was not found

    """
    if not isinstance(sums, ROOT.TList):
        raise TypeError("{0} is not of type ROOT.TList".format(sums))
    naming_pattern = "corr_this_with_{0}"
    try:
        corr_hist = asrootpy(sums.FindObject(classifier1).FindObject(naming_pattern.format(classifier2)))
    except TypeError:
        # The first classifier was not found, raising TypeError when trying to find the 2nd one
        raise ValueError("Classifier {0} was not found in given list".format(classifier1))
    if not corr_hist:
        raise ValueError("Correlation histogram {0} vs {1} not found".format(classifier1, classifier2))
    return corr_hist


def get_PNch_vs_estmult(sums, est):
    """
    Parameters
    ----------
    sums : TList
           Sums directory
    est : str
          Estimator name
    Returns
    -------
    Hist1D :
            Counter Histogram for Number of events with Nch in the estimator region
    """
    if not isinstance(sums, ROOT.TList):
        raise TypeError("{0} is not of type ROOT.TList".format(sums))
    # nasty hardcoded:
    ref_est = "EtaLt05"
    corr_hist = get_correlation_histogram(sums, est, ref_est)
    return asrootpy(corr_hist.ProjectionX(gen_random_name()))


def get_pT_distribution(results_est_dir, pids, classifier_bin_interval, normalized=False):
    """
    Parameters
    ----------
    results_est_dir : TDirectory
        Directory of a given estimator
    pids : list
        List of strings denoting requested pids
    classifier_bin_interval : tuple
        Lower and upper limit of classifier value for which the p_T distribution should be made.
        This value needs to be given as bin indices!
    normalized : Boolean
        Should the distribution be normalized to yield P(p_T)?
    Returns
    -------
    Hist1D :
        Histogram P(p_T)
    """
    mult_pt_hists = []
    for pid in pids:
        mult_pt_hists.append(getattr(results_est_dir.mult_pt, pid))
    summed_mult_pt = sum(mult_pt_hists)
    summed_mult_pt.xaxis.SetRange(*classifier_bin_interval)
    projy = asrootpy(summed_mult_pt.ProjectionY())
    projy.name = gen_random_name()
    event_counter = asrootpy(results_est_dir.event_counter)
    # Scale by the number of events in the interval;
    projy.Scale(1.0 / event_counter.Integral(*classifier_bin_interval))
    if normalized:
        projy.Scale(1.0 / projy.Integral())
    return projy


def get_mean_nMPI(sums_est_dir, classifier_bin_interval):
    """
    Get the mean nMPI of events in a given N_ch interval
    Parameters
    ----------
    sums_est_dir : TList
        List for a given estimator
    nch_low, nch_up : int
        Lower and upper limit of Nch for which <nMPI> should be calculated
    Returns
    -------
    Float :
           <nMPI>
    """
    nch_vs_nmpi = asrootpy(sums_est_dir.FindObject("corr_this_with_nMPI"))
    nch_vs_nmpi.xaxis.SetRange(*classifier_bin_interval)
    return nch_vs_nmpi.GetMean(2)


def get_graphs_particle_ratios_vs_refmult(plottingcls, pids1, pids2, scale=None, ytitle=''):
    """
    Returns list of ratios of the two pid-lists (pids1/pids2) vs refmult.
    This function depends on the correlation histograms to be present in f
    """
    ratios = []
    ref_classifier = 'EtaLt05'
    for est_dir in get_est_dirs(plottingcls.sums, plottingcls.considered_ests):
        h3d = asrootpy(est_dir.FindObject("classifier_pT_PID_{0}".format(est_dir.GetName())))
        corr_hist = get_correlation_histogram(plottingcls.sums, est_dir.GetName(), ref_classifier)
        pids1_vs_estmult = sum([get_identified_vs_mult(h3d, pdg) for pdg in pids1])
        pids2_vs_estmult = sum([get_identified_vs_mult(h3d, pdg) for pdg in pids2])

        # remap histograms using the correlation between the current estimator and the reference one
        pids1_vs_refmult = remap_x_values(pids1_vs_estmult, corr_hist)
        pids2_vs_refmult = remap_x_values(pids2_vs_estmult, corr_hist)

        # sanitize
        for g in [pids1_vs_refmult, pids2_vs_refmult]:
            remove_zero_value_points(g)
            remove_points_with_x_err_gt_1NchRef(g)
            remove_points_with_equal_x(g)
        remove_non_mutual_points(pids1_vs_refmult, pids2_vs_refmult)

        try:
            ratio = pids1_vs_refmult / pids2_vs_refmult
        except ZeroDivisionError:
            print "ZeroDivisionError in {0}".format(est_dir.GetName())
            continue
        if scale:
            ratio.Scale(scale)
        ratio.title = make_estimator_title(est_dir.GetName())
        ratios.append(ratio)
    return ratios
