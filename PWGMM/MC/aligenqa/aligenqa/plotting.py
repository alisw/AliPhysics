from pprint import pprint

from rootpy import asrootpy, log, collection
from rootpy.plotting import Hist2D
from rootpy.io import root_open

from data_extractors import \
    get_dNdeta_in_classifier_bin_interval,\
    get_identified_vs_mult,\
    get_correlation_histogram,\
    get_PNch_vs_estmult,\
    get_meanpt_vs_estmult,\
    get_pT_distribution,\
    get_mean_nMPI,\
    get_graphs_particle_ratios_vs_refmult
from utils import \
    gen_random_name,\
    get_est_dirs,\
    make_estimator_title,\
    remap_x_values,\
    remove_zero_value_points,\
    remove_points_with_equal_x,\
    remove_points_with_x_err_gt_1NchRef,\
    percentile_bin_to_binidx_bin

from .roofie import Figure, Styles

kPROTON = str(2212)
kANTIPROTON = str(-2212)
kLAMBDA = str(3122)
kANTILAMBDA = str(-3122)
kK0S = str(310)
kKPLUS = str(321)
kKMINUS = str(-321)
kPIPLUS = str(211)
kPIMINUS = str(-211)
kPI0 = str(111)
kXI = str(3312)
kANTIXI = str(-3312)
kOMEGAMINUS = str(3334)
kOMEGAPLUS = str(-3334)


class Plotting(object):
    def __init__(self, f_name, sums_dir_name, results_dir_name, percentile_bins, considered_ests):
        self.f_name = f_name
        self.sums_dir_name = sums_dir_name
        self.results_dir_name = results_dir_name
        # use the last mult bin starts at a multiplicity  x times larger than the mean in this estimator
        # self.mean_mult_cutoff_factor = 4
        self.ref_ests = ['EtaLt05', ]
        self.considered_ests = considered_ests
        self.perc_bins = percentile_bins
        # figure out the nch edges corresponding to the percentile edges, depends on P(Nch)
        self.delete_results_dir()
        self.make_results_dir()
        self.plot_event_counters()  # needed for calculations of the edges
        self.nch_edges = self._find_nch_edges_from_percentile_edges()
        pprint(self.nch_edges)
        # set the default style for all figures created from her on forward:
        Figure.style = Styles.Presentation_half

    def _io_decorator(func):
        """
        Open and close the file befor and after the execution of the decorated function.
        The purpose ist to clean up memory in this way and to force an update of the file
        before the next function calls. The wrapper adds the file, sums and results_post to `self`.
        """
        def wrapper(self, **kwargs):
            with root_open(self.f_name, 'update') as self.f:
                self.sums = self.f.MultEstimators.__getattr__(self.sums_dir_name)
                try:
                    self.results_post = self.f.MultEstimators.__getattr__(self.results_dir_name)
                except AttributeError:
                    # results dir does not exists (yet)
                    pass
                return_value = func(self, **kwargs)
                # Delete all TLists in sums since we own them and they would be left in memory otherwise
                for obj in self.sums:
                    if isinstance(obj, collection.List):
                        obj.Delete()
                self.sums.Delete()
            return return_value
        return wrapper

    @_io_decorator
    def _find_nch_edges_from_percentile_edges(self):
        nch_edges = {}
        estimators_to_be_removed = []
        for est_dir in get_est_dirs(self.results_post, self.considered_ests):
            event_counter = est_dir.event_counter
            try:
                nch_edges[est_dir.GetName()] = [percentile_bin_to_binidx_bin(perc_bin, event_counter)
                                                for perc_bin in self.perc_bins[est_dir.GetName()]]
            except ValueError, e:
                print "Error occured for classifier " + est_dir.GetName()
                print e
                print self.perc_bins[est_dir.GetName()]
                print "You can change the percentile bins in the beginning of this script"
                print "For the following, this estimator is removed"
                estimators_to_be_removed.append(est_dir.GetName())
        print "Bin edges for given percentile bins"
        print nch_edges
        for est in estimators_to_be_removed:
            del self.perc_bins[est]
            del self.considered_ests[self.considered_ests.index(est)]
        return nch_edges

    @_io_decorator
    def delete_results_dir(self):
        # delete old result directory
        self.f.rm('MultEstimators/' + self.results_dir_name)
        self.f.Write()

    @_io_decorator
    def make_results_dir(self):
        self.f.mkdir('MultEstimators/' + self.results_dir_name, recurse=True)
        for est_dir in get_est_dirs(self.sums, self.considered_ests):
            try:
                resdir = self.f.MultEstimators.__getattr__(self.results_dir_name).mkdir(est_dir.GetName())
                resdir.Write()
            except:
                pass

    @_io_decorator
    def plot_particle_ratios_vs_estmult(self, pids1, pids2, scale=None, ytitle=''):
        ratio_vs_estmult_dir = (self.results_post.GetPath().split(":")[1]  # file.root:/internal/root/path
                                + '/pid_ratios_vs_estmult')
        fig = Figure()
        if not ytitle:
            fig.ytitle = ", ".join(pids1) + " / " + ", ".join(pids2)
        else:
            fig.ytitle = ytitle

        for est_dir in get_est_dirs(self.sums, self.considered_ests):
            h3d = asrootpy(est_dir.FindObject("fNch_pT_pid"))
            pids1hists = [get_identified_vs_mult(h3d, pdg) for pdg in pids1]
            pids2hists = [get_identified_vs_mult(h3d, pdg) for pdg in pids2]

            pids1_px = sum(pids1hists)
            pids2_px = sum(pids2hists)
            ratio1d = pids1_px / pids2_px

            fig.xtitle = "N_{ch}|_{" + make_estimator_title(est_dir.GetName()) + "}"

            if scale:
                ratio1d.Scale(scale)
                fig.add_plottable(ratio1d, legend_title=make_estimator_title(est_dir.GetName()))
                name = "_".join(pids1) + "_div_" + "_".join(pids2)
                fig.save_to_root_file(self.f, name, ratio_vs_estmult_dir)

    @_io_decorator
    def plot_event_counters(self):
        log.info("Creating event counters")
        for est_dir in get_est_dirs(self.sums, self.considered_ests):
            results_est_dir = self.results_post.__getattr__(est_dir.GetName())
            # Nasty, but just use a reference estimator here...
            corr = get_correlation_histogram(self.sums, est_dir.GetName(), "EtaLt05")
            counter = asrootpy(corr.ProjectionX())
            counter.name = "event_counter"
            path = results_est_dir.GetPath().split(":")[1]  # file.root:/internal/root/path
            self.f.cd(path)
            results_est_dir.WriteTObject(counter)

    @_io_decorator
    def plot_dNdetas(self, ratio_to_mb):
        # Loop over all estimators in the Sums list:
        log.info("Creating dN/deta bin in multiplicity")
        figs = []
        for est_dir in get_est_dirs(self.sums, self.considered_ests):
            # does this estimator have several multiplicity bins?
            # Q2, for example only works with pythia and makes no sense to plot
            # on Dipsy as it would only be the MB line
            if len(self.nch_edges[est_dir.GetName()]) == 1:
                continue
            results_est_dir = self.results_post.Get(est_dir.GetName())
            event_counter = asrootpy(results_est_dir.Get("event_counter"))

            fig = Figure()
            fig.plot.palette = 'colorblind'
            fig.xtitle = '#eta'
            fig.ytitle = 'Ratio of dN_{ch}/d#eta over MB result' if ratio_to_mb else '1/N #times dN_{ch}/d#eta'
            fig.legend.title = make_estimator_title(est_dir.GetName())
            fig.plot.ymin = 0
            dNdeta_mb = get_dNdeta_in_classifier_bin_interval(est_dir, event_counter,
                                                              [1, event_counter.GetXaxis().GetNbins()])

            for cls_bin, perc_bin in zip(self.nch_edges[est_dir.GetName()], self.perc_bins[est_dir.GetName()]):
                title = "{0}%-{1}%".format(perc_bin[1] * 100, perc_bin[0] * 100)
                dNdeta_in_interval = get_dNdeta_in_classifier_bin_interval(est_dir, event_counter, cls_bin)
                if ratio_to_mb:
                    fig.add_plottable(dNdeta_in_interval / dNdeta_mb, legend_title=title)
                else:
                    fig.add_plottable(dNdeta_in_interval, legend_title=title)
            # add MB as well, if it is not the ratio plots we are making
            if not ratio_to_mb:
                title = "MB"
                fig.add_plottable(dNdeta_mb, legend_title=title)
            path = results_est_dir.GetPath().split(":")[1]  # file.root:/internal/root/path
            if ratio_to_mb:
                fig.save_to_root_file(self.f, "dNdeta_MB_ratio_summary", path=path)
            else:
                fig.save_to_root_file(self.f, "dNdeta_summary", path=path)
            figs.append(fig)
        return figs

    @_io_decorator
    def plot_pt_distribution_ratios(self):
        # create particle ratio vs pT plots
        log.info("Computing histograms vs pt")
        results_path = self.results_post.GetPath().split(":")[1]  # file.root:/internal/root/path
        # Loop over all estimators in the Sums list:
        figs = []

        def get_new_figure():
            fig = Figure()
            fig.xtitle = 'p_{T} (GeV)'
            fig.plot.ymin = 0
            fig.plot.xmax = 10
            fig.plot.palette = 'colorblind'
            # fig.plot.palette_ncolors = len(nch_edges) - 1
            fig.legend.position = 'br'
            return fig

        for est_dir in get_est_dirs(self.results_post, self.considered_ests):
            dirname = '{0}/{1}/pid_ratios/'.format(results_path, est_dir.GetName())

            mult_binned_pt_dists = {}
            mult_binned_pt_dists['proton'] = [
                get_pT_distribution(est_dir, [kANTIPROTON, kPROTON], classifier_bin_interval)
                for classifier_bin_interval in self.nch_edges[est_dir.GetName()]
            ]
            mult_binned_pt_dists['pi_ch'] = [
                get_pT_distribution(est_dir, [kPIMINUS, kPIPLUS], classifier_bin_interval)
                for classifier_bin_interval in self.nch_edges[est_dir.GetName()]
            ]
            mult_binned_pt_dists['xi'] = [
                get_pT_distribution(est_dir, [kANTIXI, kXI], classifier_bin_interval)
                for classifier_bin_interval in self.nch_edges[est_dir.GetName()]
            ]
            mult_binned_pt_dists['omega'] = [
                get_pT_distribution(est_dir, [kOMEGAMINUS, kOMEGAPLUS], classifier_bin_interval)
                for classifier_bin_interval in self.nch_edges[est_dir.GetName()]
            ]
            mult_binned_pt_dists['lambda'] = [
                get_pT_distribution(est_dir, [kANTILAMBDA, kLAMBDA], classifier_bin_interval)
                for classifier_bin_interval in self.nch_edges[est_dir.GetName()]
            ]
            mult_binned_pt_dists['k0s'] = [
                get_pT_distribution(est_dir, [kK0S], classifier_bin_interval)
                for classifier_bin_interval in self.nch_edges[est_dir.GetName()]
            ]
            mult_binned_pt_dists['k_ch'] = [
                get_pT_distribution(est_dir, [kKPLUS, kKMINUS], classifier_bin_interval)
                for classifier_bin_interval in self.nch_edges[est_dir.GetName()]
            ]
            mult_binned_pt_dists['pi0'] = [
                get_pT_distribution(est_dir, [kPI0], classifier_bin_interval)
                for classifier_bin_interval in self.nch_edges[est_dir.GetName()]
            ]
            perc_titles = ["{0}%-{1}%".format(perc_bin[1] * 100, perc_bin[0] * 100)
                           for perc_bin in self.perc_bins[est_dir.GetName()]]

            fig = get_new_figure()
            name = "proton_over_pich__vs__pt"
            fig.ytitle = "(p+#bar{p})/#pi^{+-}"
            fig.plot.ymax = .3
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['proton'], mult_binned_pt_dists['pi_ch'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            fig = get_new_figure()
            name = "Xi_over_pich__vs__pt"
            fig.plot.ymax = .06
            fig.legend.position = 'tl'
            fig.ytitle = "#Xi/#pi^{+-}"
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['xi'], mult_binned_pt_dists['pi_ch'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            fig = get_new_figure()
            name = "OmegaCh_over_pich__vs__pt"
            fig.plot.ymax = .005
            fig.legend.position = 'tl'
            fig.ytitle = "#Omega_{ch}/#pi^{+-} "
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['omega'], mult_binned_pt_dists['pi_ch'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            # Ratios to pi0
            fig = get_new_figure()
            name = "pich_over_pi0__vs__pt"
            fig.plot.ymax = 2.5
            fig.legend.position = 'bl'
            fig.ytitle = "#pi^{+-}/#pi^{0}"
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['pi_ch'], mult_binned_pt_dists['pi0'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            fig = get_new_figure()
            name = "proton_over_pi0__vs__pt"
            fig.plot.ymax = 1
            fig.legend.position = 'tr'
            fig.ytitle = "p/#pi^{0}"
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['proton'], mult_binned_pt_dists['pi0'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            fig = get_new_figure()
            name = "K0S_over_pi0__vs__pt"
            fig.plot.ymax = 1.4
            fig.legend.position = 'tl'
            fig.ytitle = "K^{0}_{S}/#pi^{0}"
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['k0s'], mult_binned_pt_dists['pi0'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            fig = get_new_figure()
            name = "Lambda_over_pi0__vs__pt"
            fig.plot.ymax = .9
            fig.legend.position = 'tl'
            fig.ytitle = "#Lambda/#pi^{0}"
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['lambda'], mult_binned_pt_dists['pi0'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            fig = get_new_figure()
            name = "Xi_over_pi0__vs__pt"
            fig.plot.ymax = .08
            fig.legend.position = 'tl'
            fig.ytitle = "#Xi/#pi^{0}"
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['xi'], mult_binned_pt_dists['pi0'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            fig = get_new_figure()
            name = "OmegaCh_over_pi0__vs__pt"
            fig.plot.ymax = .005
            fig.legend.position = 'tl'
            fig.ytitle = "#Omega_{ch}/#pi^{0}"
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['omega'], mult_binned_pt_dists['pi0'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            # Ratios to K0S
            fig = get_new_figure()
            name = "proton_over_K0S__vs__pt"
            fig.plot.ymax = 2.6
            fig.legend.position = 'tr'
            fig.ytitle = "p/K^{0}_{S}"
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['proton'], mult_binned_pt_dists['k0s'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            fig = get_new_figure()
            name = "Lambda_over_K0S__vs__pt"
            fig.plot.ymax = 1
            fig.legend.position = 'bl'
            fig.ytitle = "#Lambda/K^{0}_{S}"
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['lambda'], mult_binned_pt_dists['k0s'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            fig = get_new_figure()
            name = "Xi_over_K0S__vs__pt"
            fig.plot.ymax = .2
            fig.legend.position = 'tl'
            fig.ytitle = "#Xi/K^{0}_{S}"
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['xi'], mult_binned_pt_dists['k0s'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            fig = get_new_figure()
            name = "OmegaCh_over_K0S__vs__pt"
            fig.plot.ymax = .012
            fig.legend.position = 'tl'
            fig.ytitle = "#Omega_{ch}/K^{0}_{S}"
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['omega'], mult_binned_pt_dists['k0s'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

            fig = get_new_figure()
            name = "Kaon_over_pich__vs__pt"
            fig.plot.ymax = 1
            fig.legend.position = 'tl'
            fig.ytitle = "(K^{+} + K^{-}) / (#pi^{+} +#pi^{-})"
            fig.legend.title = make_estimator_title(est_dir.GetName())
            [
                fig.add_plottable(h1 / h2, legend_title=title)
                for h1, h2, title in zip(mult_binned_pt_dists['k_ch'], mult_binned_pt_dists['pi_ch'], perc_titles)
            ]
            fig.save_to_root_file(self.f, name, dirname)
            figs.append(fig)

        return figs

    @_io_decorator
    def plot_PNch_summary(self):
        log.info("Creating P(Nch) summary plot")
        summary_fig = Figure()
        summary_fig.xtitle = "N_{ch}^{est}"
        summary_fig.ytitle = "P(N_{ch}^{est})"
        summary_fig.legend.position = 'tr'
        summary_fig.plot.logy = True

        for est_dir in get_est_dirs(self.sums, self.considered_ests):
            est_name = est_dir.GetName()
            h_tmp = get_PNch_vs_estmult(self.sums, est_name)
            if h_tmp.Integral() > 0:
                h_tmp.Scale(1.0 / h_tmp.Integral())
                summary_fig.add_plottable(h_tmp, make_estimator_title(est_name))

        path = self.results_post.GetPath().split(":")[1]  # file.root:/internal/root/path
        summary_fig.save_to_root_file(self.f, "PNch_summary", path=path)
        # list as return type is expected for making the pdf
        return [summary_fig]

    @_io_decorator
    def plot_PNch(self):
        log.info("Creating P(Nch_est) and P(Nch_refest) histograms")
        # mult_bin_size = 10
        figs = []
        for ref_est_name in self.ref_ests:
            for res_est_dir in get_est_dirs(self.results_post, self.considered_ests):
                est_name = res_est_dir.GetName()
                # Figure properties:
                fig_vs_estmult = Figure()
                fig_vs_refmult = Figure()
                fig_vs_estmult.plot.logy = True
                fig_vs_refmult.plot.logy = True

                fig_vs_estmult.plot.palette = 'colorblind'
                fig_vs_refmult.plot.palette = 'colorblind'

                fig_vs_estmult.legend.position = 'tr'
                fig_vs_refmult.legend.position = 'tr'

                fig_vs_estmult.xtitle = "N_{{ch}}^{{{0}}}".format(est_name)
                fig_vs_refmult.xtitle = "N_{{ch}}^{{{0}}}".format(ref_est_name)

                fig_vs_estmult.ytitle = "P(N_{{ch}}^{{{0}}})".format(est_name)
                fig_vs_refmult.ytitle = "P(N_{{ch}}^{{{0}}})".format(ref_est_name)

                corr_hist = get_correlation_histogram(self.sums, est_name, ref_est_name)

                # logic when dealing with fixed bins given in Nch:
                # ------------------------------------------------
                # mean_nch_est = corr_hist.GetMean(1)  # mean of x axis
                # nch_max = corr_hist.xaxis.GetNbins()
                # nch_cutoff = mean_nch_est * mean_mult_cutoff_factor
                # nch_bins = [(low, low + mult_bin_size) for low in range(0, int(nch_cutoff), mult_bin_size)]
                # # a large last bin covering the rest:
                # nch_bins += [(nch_bins[-1][2], nch_max)]
                # legend_tmpl = "{} < N_{ch} < {}"
                # logic when dealing with percentile bins:
                # ----------------------------------------
                # event_counter_est = asrootpy(getattr(res_est_dir, "event_counter"))

                legend_tmpl = "{0}% - {1}%"
                fig_vs_estmult.legend.title = "Selected in {0}".format(make_estimator_title(ref_est_name))
                fig_vs_refmult.legend.title = "Selected in {0}".format(make_estimator_title(est_name))
                # WARNING: the following needs tweeking when going back to fixed N_ch bins!
                for nch_bin, perc_bin in zip(self.nch_edges[ref_est_name], self.perc_bins[ref_est_name]):
                    # vs est_mult:
                    corr_hist.xaxis.SetRange(0, 0)  # reset x axis
                    corr_hist.yaxis.SetRange(nch_bin[0], nch_bin[1])
                    h_vs_est = asrootpy(corr_hist.ProjectionX(gen_random_name()))
                    if h_vs_est.Integral() > 0:
                        h_vs_est.Scale(1.0 / h_vs_est.Integral())
                        fig_vs_estmult.add_plottable(h_vs_est, legend_tmpl.format(perc_bin[1] * 100, perc_bin[0] * 100))
                    else:
                        log.info("No charged particles in {0}*100 percentile bin of estimator {1}. This should not happen".
                                 format(perc_bin, ref_est_name))
                for nch_bin, perc_bin in zip(self.nch_edges[est_name], self.perc_bins[est_name]):
                    # vs ref_mult:
                    corr_hist.yaxis.SetRange(0, 0)  # reset y axis
                    corr_hist.xaxis.SetRange(*nch_bin)
                    h_vs_ref = asrootpy(corr_hist.ProjectionY(gen_random_name()))
                    if h_vs_ref.Integral() > 0:
                        h_vs_ref.Scale(1.0 / h_vs_ref.Integral())
                        fig_vs_refmult.add_plottable(h_vs_ref, legend_tmpl.format(perc_bin[1] * 100, perc_bin[0] * 100))
                    else:
                        log.info(
                            "No charged particles in {0}*100 percentile bin of estimator {1}. This should not happen".
                            format(perc_bin, est_name))

                path = res_est_dir.GetPath().split(":")[1]
                # vs est_mult
                fig_vs_estmult.save_to_root_file(self.f, "PNchEst_binned_in_Nch{0}".format(ref_est_name), path)
                # vs est_mult
                fig_vs_refmult.save_to_root_file(self.f, "PNch{0}_binned_in_NchEst".format(ref_est_name), path)
                figs.append(fig_vs_estmult)
                figs.append(fig_vs_refmult)
        return figs

    @_io_decorator
    def plot_mult_vs_pt(self):
        log.info("Makeing 2D  pt plots for each particle kind")
        for est_dir in get_est_dirs(self.sums, self.considered_ests):
            path = (self.results_post.GetPath().split(":")[1]  # file.root:/internal/root/path
                    + "/" + est_dir.GetName()
                    + "/mult_pt")
            try:
                self.f.mkdir(path, recurse=True)
            except ValueError:
                pass
            self.f.cd(path)

            h3d = asrootpy(est_dir.FindObject('classifier_pT_PID_{0}'.format(est_dir.GetName())))
            # loop through all particle kinds:
            nPIDs = h3d.zaxis.GetNbins()
            for ibin in range(1, nPIDs + 1):
                h3d.zaxis.SetRange(ibin, ibin)
                mult_pt = asrootpy(h3d.Project3D("yx"))
                mult_pt.name = h3d.zaxis.GetBinLabel(ibin)
                mult_pt.Write()

    @_io_decorator
    def plot_correlation(self):
        # Make correlations between estimators
        log.info("Correlating N_ch of each estimator")
        corr_dir = self.results_post.GetPath().split(":")[1] + '/correlations'
        try:
            self.f.mkdir(corr_dir, recurse=True)
        except:
            pass
        # Take ntuple from the first estimator and then add friends to this one
        nt0 = self.sums[0].FindObject("fEventTuple")
        nt0.SetAlias(self.sums[0].GetName(), "fEventTuple")

        # build ntuple
        for est_dir in self.sums[1:]:
            nt0.AddFriend(est_dir.FindObject("fEventTuple"), est_dir.GetName())
        for ref_est in self.considered_ests:
            for est_dir in self.sums:
                log.info("Correlating {0} with {1}".format(ref_est, est_dir.GetName()))
                corr_hist = Hist2D(400, 0, 400,
                                   400, 0, 400,
                                   name="corr_hist_{0}_vs_{1}".format(ref_est, est_dir.GetName()))
                # Lables are deliberatly swaped, see Projection below!
                corr_hist.title = ("Correlation N_{{ch}} in {0} and {1};N_{{ch}} {1};N_{{ch}} {0}"
                                   .format(ref_est, est_dir.GetName()))

                # this projects onto y:x, to make coding more adventurous
                nt0.Project(corr_hist.name, "{0}.nch:{1}.nch".format(ref_est, est_dir.GetName()),
                            "ev_weight")
                corr_hist.drawstyle = 'colz'
                self.f.cd(corr_dir)
                corr_hist.write()

    @_io_decorator
    def plot_pid_ratio_vs_refmult(self):
        log.info("Creating plots vs refmult")
        ratios_dir = self.results_post.GetPath().split(":")[1] + '/pid_ratios_vs_refmult'

        def get_new_figure():
            fig = Figure()
            fig.plot.ncolors = len(self.considered_ests)
            fig.xtitle = "N_{ch}|_{" + make_estimator_title('EtaLt05') + "}"
            fig.plot.xmin = 0
            fig.plot.xmax = 60
            return fig

        figs = []
        # Proton / pi_ch
        fig = get_new_figure()
        pids1, pids2 = ['-2212', '2212'], ['-211', '211']
        fig.ytitle = "p/#pi^{+-}"
        fig.plot.ymin, fig.plot.ymax = 0.04, 0.13
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2, )
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # K / pi_ch
        fig = get_new_figure()
        pids1, pids2 = ['310', '321', '-321'], ['-211', '211']
        fig.ytitle = "K^{*}/#pi^{+-}"
        fig.plot.ymin, fig.plot.ymax = 0.09, 0.30
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # Lambda / pi_ch
        fig = get_new_figure()
        pids1, pids2 = ['3122'], ['-211', '211']
        fig.ytitle = "#Lambda / #pi^{+-}"
        fig.plot.ymin, fig.plot.ymax = 0.005, 0.035
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # Xi / pi_ch
        fig = get_new_figure()
        pids1, pids2 = ['3312'], ['-211', '211']
        fig.ytitle = "#Xi / #pi^{+-}"
        fig.plot.ymin, fig.plot.ymax = 0.0004, 0.003
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # Omega / pi_ch
        fig = get_new_figure()
        pids1, pids2 = ['3334', '-3334'], ['-211', '211']
        fig.ytitle = "#Omega / #pi^{+-}"
        fig.plot.ymin, fig.plot.ymax = 0.00001, 0.0005
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # pi_ch/pi0
        fig = get_new_figure()
        pids1, pids2 = ['-211', '211'], ['111']
        fig.ytitle = "#pi^{+-}/#pi^{0}"
        fig.plot.ymin, fig.plot.ymax = 1.5, 2.2
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # proton / pi0
        fig = get_new_figure()
        pids1, pids2 = ['-2212', '2212'], ['111']
        fig.ytitle = "p/#pi^{0}"
        fig.plot.ymin, fig.plot.ymax = 0.09, 0.30
        fig.legend.position = 'tl'
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # K / pi0
        fig = get_new_figure()
        pids1, pids2 = ['310', '321', '-321'], ['111']
        fig.ytitle = "K^{*}/#pi^{0}"
        fig.plot.ymin, fig.plot.ymax = 0.15, 0.50
        fig.legend.position = 'tl'
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # Lambda / pi0
        fig = get_new_figure()
        pids1, pids2 = ['3122'], ['111']
        fig.ytitle = "#Lambda/#pi^{0}"
        fig.plot.ymin, fig.plot.ymax = 0.014, 0.045
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # Xi / pi0
        fig = get_new_figure()
        pids1, pids2 = ['3312'], ['111']
        fig.ytitle = "#Xi/#pi^{0}"
        fig.plot.ymin, fig.plot.ymax = 0.0010, 0.005
        fig.legend.position = 'tl'
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # Omega / pi0
        fig = get_new_figure()
        pids1, pids2 = ['3334', '-3334'], ['111']
        fig.ytitle = "#Omega/#pi^{0}"
        fig.legend.position = 'tl'
        fig.plot.ymin, fig.plot.ymax = 0.00002, 0.0008
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # K_ch / K0_S
        fig = get_new_figure()
        pids1, pids2 = ['321', '-321'], ['310']
        fig.ytitle = "(K^{+}+K^{-}) / (2#timesK^{0}_{S})"
        fig.plot.ymin, fig.plot.ymax = 0.4, 1.5
        fig.legend.position = 'tl'
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2, scale=.5)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # K0_S / Lambda
        fig = get_new_figure()
        pids1, pids2 = ['310'], ['-3122', '3122']
        fig.ytitle = "K^{0}_{S} / #Lambda"
        fig.plot.ymin, fig.plot.ymax = 1.3, 3.7
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        # K0_S / Xi
        fig = get_new_figure()
        pids1, pids2 = ['310'], ['3312']
        fig.ytitle = "K^{0}_{S} / #Xi"
        fig.plot.ymin, fig.plot.ymax = 15, 80
        graphs = get_graphs_particle_ratios_vs_refmult(self, pids1, pids2)
        [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
        name = "_".join(pids1) + "_div_" + "_".join(pids2)
        fig.save_to_root_file(self.f, name, ratios_dir)
        figs.append(fig)

        return figs

        # ######################################################################################
        # # vs Est mult
        # _plot_particle_ratios_vs_estmult(self, ['321', '-321'], ['310'],
        #                                  scale=.5, fig.ytitle = "(K^{+} + K^{-}) / (2*K_{S}^{0})")

    @_io_decorator
    def plot_meanpt_vs_ref_mult_for_pids(self):
        log.info("Creating mean pT plots")
        figs = []
        for sums_est_dir, res_est_dir in zip(get_est_dirs(self.sums, self.considered_ests),
                                             get_est_dirs(self.results_post, self.considered_ests)):
            if sums_est_dir.GetName() != res_est_dir.GetName():
                raise IndexError("Order of estimator dirs is different in sums and results_post")
            res_dir_str = res_est_dir.GetPath().split(":")[1]
            corr_hist = get_correlation_histogram(self.sums, sums_est_dir.GetName(), "EtaLt05")
            # Get the <pT> per classifier bin; then, re-map the classifier value to the reference classifier (eg EtaLt05)
            # This might not make a lot of sense, actually. Maybe it would be much more telling if I were to
            # put the percentile bins on the x-axis? As in the highest 1% of that classifier has a <pT> of ...
            graphs = []
            graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kPI0, kPIMINUS, kPIPLUS]), corr_hist))
            graphs[-1].title = "#pi"
            graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kKMINUS, kKPLUS]), corr_hist))
            graphs[-1].title = "K^{#pm}"
            graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kPROTON, kANTIPROTON]), corr_hist))
            graphs[-1].title = "p"
            graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kK0S]), corr_hist))
            graphs[-1].title = "K^{0}_{S}"
            graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kLAMBDA, kANTILAMBDA]), corr_hist))
            graphs[-1].title = "#Lambda"
            graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kXI, kANTIXI]), corr_hist))
            graphs[-1].title = "#Xi"
            graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kOMEGAMINUS, kOMEGAPLUS]), corr_hist))
            graphs[-1].title = "#Omega"
            # sanitize graphs:
            for g in graphs:
                remove_zero_value_points(g)
                remove_points_with_x_err_gt_1NchRef(g)
                remove_points_with_equal_x(g)

            fig = Figure()
            fig.plot.palette = 'root'
            fig.plot.ncolors = 7
            fig.plot.xmin = 0
            fig.plot.xmax = 40
            fig.plot.ymin = 0.3
            fig.plot.ymax = 2.1
            fig.ytitle = "<p_{T}>"
            fig.xtitle = "N_{ch}|_{|#eta|<0.5}"
            fig.legend.title = make_estimator_title(sums_est_dir.GetName())
            [fig.add_plottable(g, g.title) for g in graphs]
            fig.save_to_root_file(self.f, "mean_pt", res_dir_str)
            figs.append(fig)
        return figs

    # def _plot_event_counter_with_shaded_perc_areas(f, results_post):
    #     log.info("Broken: Root sucks! Creating shaded event counter with percentile regions")
    #     return
    #     for est_dir in get_est_dirs(results_post):
    #         event_counter = asrootpy(getattr(est_dir, "event_counter"))
    #         nch_edges = get_Nch_edges_for_percentile_edges(perc_edges, event_counter)
    #         c = Canvas(name="event_counter_with_perc")
    #         leg = Legend(len(nch_edges) - 1)
    #         copies = []
    #         colors = get_color_generator(ncolors=10)
    #         # Draw the hist once
    #         event_counter.Draw()
    #         for nch_low, nch_up in zip(nch_edges[:-1], nch_edges[1:]):
    #             copies.append(event_counter.Clone(gen_random_name()))
    #             copies[-1].xaxis.SetRangeUser(nch_low, nch_up)
    #             copies[-1].SetFillStyle(1001)
    #             copies[-1].color = next(colors)
    #             copies[-1].xaxis.title = "N_{ch}"
    #             copies[-1].yaxis.title = "counts"
    #             leg.AddEntry(copies[-1], "{}-{}%".format(str(nch_low), str(nch_up)))
    #             copies[-1].Draw('sameHist')
    #             break
    #         leg.Draw()
    #         est_dir.cd()
    #         c.Write()

    @_io_decorator
    def plot_dNdpT(self, pid_selection):
        """
        Plot dNdpT particles in pid_selection

        Parameters
        ----------
        pid_selection : str
            Either all charged particles ('ch') or 'pi', 'K' or 'p'
        """
        log.info("1/N_evts  dN_ch/dpT plots")
        figs = []
        for sums_est_dir, res_est_dir in zip(get_est_dirs(self.sums, self.considered_ests),
                                             get_est_dirs(self.results_post, self.considered_ests)):
            if sums_est_dir.GetName() != res_est_dir.GetName():
                raise IndexError("Order of estimator dirs is different in sums and results_post")
            res_dir_str = res_est_dir.GetPath().split(":")[1]
            fig = Figure()
            fig.plot.palette = 'colorblind'
            # fig.plot.ncolors = 5
            fig.legend.position = 'tr'
            fig.ytitle = "1/N_{evts} dN/dp_{T} (" + make_estimator_title(sums_est_dir.GetName()) + ")"
            fig.xtitle = "p_{T} (GeV)"
            fig.plot.logy = True
            hists = []
            if pid_selection == 'ch':
                fig.legend.title = "#pi^{#pm}, K^{#pm}, p, #Lambda, #Xi, #Omega"
                pid_numbers = [kPIMINUS, kPIPLUS, kKMINUS, kKPLUS, kPROTON, kANTIPROTON,
                               kLAMBDA, kANTILAMBDA, kXI, kANTIXI, kOMEGAMINUS, kOMEGAPLUS]
            if pid_selection == 'pi':
                fig.legend.title = "#pi^{#pm}"
                pid_numbers = [kPIMINUS, kPIPLUS]
            if pid_selection == 'K':
                fig.legend.title = "K^{#pm}"
                pid_numbers = [kKMINUS, kKPLUS]
            if pid_selection == 'p':
                fig.legend.title = "p, #bar{p}"
                pid_numbers = [kPROTON, kANTIPROTON]

            for perc_bin, classifier_bin in zip(self.perc_bins[sums_est_dir.GetName()], self.nch_edges[sums_est_dir.GetName()]):
                hists.append(get_pT_distribution(res_est_dir, pid_numbers, classifier_bin, normalized=False))
                hists[-1].title = "{0}%-{1}%".format(perc_bin[1] * 100, perc_bin[0] * 100)

            # add MB last to be consistent with colors in other plots; the very first and very last bin we look at
            classifier_bin_mb = (self.nch_edges[sums_est_dir.GetName()][0][0], self.nch_edges[sums_est_dir.GetName()][-1][-1])
            hists.append(get_pT_distribution(res_est_dir, pid_numbers, classifier_bin_mb, normalized=False))
            hists[-1].title = "MB"

            # scale by bin width
            [h.Scale(1, "width") for h in hists]

            [fig.add_plottable(p, p.title) for p in hists]
            fig.save_to_root_file(self.f, "dN{0}dpT".format(pid_selection), res_dir_str)
            figs.append(fig)
        return figs

    @_io_decorator
    def plot_pT_HM_div_pt_MB(self, scale_nMPI):
        log.info("Plot dN_{HM}/dpT / dN_{MB}/dpT ratios scaled with nMPI")
        figs = []
        for sums_est_dir, res_est_dir in zip(get_est_dirs(self.sums, self.considered_ests),
                                             get_est_dirs(self.results_post, self.considered_ests)):
            if sums_est_dir.GetName() != res_est_dir.GetName():
                raise IndexError("Order of estimator dirs is different in sums and results_post")
            res_dir_str = res_est_dir.GetPath().split(":")[1]
            fig = Figure()
            fig.plot.palette = 'root'
            fig.plot.ncolors = 7
            fig.xtitle = "p_{T} (GeV)"
            fig.legend.title = make_estimator_title(sums_est_dir.GetName())
            if scale_nMPI:
                fig.ytitle = ("#left[ #frac{dN^{HM}}{dp_{T}} / #frac{dN^{MB}}{dp_{T}} #right] "
                              "#times #left[ #frac{<N_{MPI}^{MB}>}{<N_{MPI}^{HM}>} #right]")
            else:
                fig.ytitle = "#frac{dN^{HM}}{dp_{T}} / #frac{dN^{MB}}{dp_{T}}"

            charged_particles = [kPIMINUS, kPIPLUS, kKMINUS, kKPLUS, kPROTON, kANTIPROTON,
                                 kLAMBDA, kANTILAMBDA, kXI, kANTIXI, kOMEGAMINUS, kOMEGAPLUS]

            # get the MB distribution which will be used to devide the nch-binned distributions
            classifier_bin_mb = (self.nch_edges[sums_est_dir.GetName()][0][0],
                                 self.nch_edges[sums_est_dir.GetName()][-1][-1])
            pt_dist_mb = get_pT_distribution(res_est_dir, charged_particles, classifier_bin_mb, normalized=False)
            mean_nmpi_mb = get_mean_nMPI(sums_est_dir, classifier_bin_mb)

            for perc_bin, classifier_bin in zip(self.perc_bins[sums_est_dir.GetName()],
                                                self.nch_edges[sums_est_dir.GetName()]):
                # get the pt distribution in this Nch interval
                pt_dist_in_interval = get_pT_distribution(res_est_dir, charged_particles,
                                                          classifier_bin, normalized=False)
                title = "{0}%-{1}%".format(perc_bin[1] * 100, perc_bin[0] * 100)
                if scale_nMPI:
                    mean_nmpi_hm = get_mean_nMPI(sums_est_dir, classifier_bin)
                    fig.add_plottable((pt_dist_in_interval / pt_dist_mb) * (mean_nmpi_mb / mean_nmpi_hm), title)
                    name = "pt_hm_div_pt_mb_scaled_nMPI"
                else:
                    fig.add_plottable((pt_dist_in_interval / pt_dist_mb), title)
                    name = "pt_hm_div_pt_mb"
            fig.save_to_root_file(self.f, name, res_dir_str)
            figs.append(fig)
        return figs

    @_io_decorator
    def plot_nMPI_vs_Nch(self):
        log.info("Creating nMPI(Nch) summary plot")
        summary_fig = Figure()
        summary_fig.xtitle = "N_{ch}^{est}"
        summary_fig.ytitle = "<N_{MPI}>"
        summary_fig.plot.palette = 'root'
        summary_fig.legend.position = 'br'
        summary_fig.plot.logy = True
        summary_fig.plot.ymin = 1

        for est_dir in get_est_dirs(self.sums, self.considered_ests):
            h_tmp = asrootpy(get_correlation_histogram(self.sums, est_dir.GetName(), "nMPI").ProfileX())
            summary_fig.add_plottable(h_tmp, make_estimator_title(est_dir.GetName()))

        path = self.results_post.GetPath().split(":")[1]  # file.root:/internal/root/path
        summary_fig.save_to_root_file(self.f, "nMPI_summary", path=path)
        return [summary_fig]
