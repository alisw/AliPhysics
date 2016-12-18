import ROOT
import settings

from rootpy.io import root_open

from hmtfmc import roofie, utils


def summarize(args):
    """
    Create summary slides for all triggers of a given
    generator. Remember to use quotes if the generator name contains
    spaces and be careful with special characters.
    """
    ROOT.gROOT.SetBatch(True)

    for global_trigger in settings.considered_triggers:
        with root_open(args.input_file, "read") as f:
            results_dir_name = "results_post" + global_trigger

            if args.gen_name:
                gen_name = args.gen_name
            else:
                gen_name = utils.get_generator_name_from_filename(args.input_file)
            latexdoc = roofie.Beamerdoc(author="HMTF (Christian Bourjau)",
                                        title=r"{0} {1}".format(gen_name, global_trigger))

            # Fish the plots we want out of the .root file and if necessary adjust some visual settings
            sec = latexdoc.add_section(r"Highlights")

            c = f.MultEstimators.__getattr__(results_dir_name).__getattr__("EtaLt15").__getattr__('dNdeta_summary')
            sec.add_figure(c)
            c = f.MultEstimators.__getattr__(results_dir_name).__getattr__("V0M").__getattr__('dNdeta_summary')
            sec.add_figure(c)
            c = f.MultEstimators.__getattr__(results_dir_name).__getattr__("EtaLt15").__getattr__('dNdeta_MB_ratio_summary')
            sec.add_figure(c)
            c = f.MultEstimators.__getattr__(results_dir_name).__getattr__("V0M").__getattr__('dNdeta_MB_ratio_summary')
            sec.add_figure(c)

            # this plot needs extra massaging, we import it to a roofie figure
            c = f.MultEstimators.__getattr__(results_dir_name).__getattr__('PNch_summary')
            fig_all = roofie.Figure()
            fig_sub = roofie.Figure()
            fig_all.import_plottables_from_canvas(c)
            # filter out an reduced set of estimators
            for p in fig_all._plottables:
                if p['legend_title'] in [r'|#eta|#leq1.5', 'ZDC', 'V0M']:
                    fig_sub.add_plottable(p['p'], p['legend_title'])
            # configure the plot
            fig_sub.xtitle = fig_all.xtitle
            fig_sub.ytitle = fig_all.ytitle
            fig_sub.plot.logy = True
            fig_sub.legend.position = 'tr'
            sec.add_figure(fig_sub)

            c = f.MultEstimators.__getattr__(results_dir_name).EtaLt15.dNchdpT
            sec.add_figure(c)
            c = f.MultEstimators.__getattr__(results_dir_name).EtaLt15.dNpdpT
            sec.add_figure(c)
            c = f.MultEstimators.__getattr__(results_dir_name).EtaLt15.dNpidpT
            sec.add_figure(c)
            c = f.MultEstimators.__getattr__(results_dir_name).EtaLt15.dNKdpT
            sec.add_figure(c)

            c = f.MultEstimators.__getattr__(results_dir_name).V0M.dNchdpT
            sec.add_figure(c)
            c = f.MultEstimators.__getattr__(results_dir_name).V0M.dNpdpT
            sec.add_figure(c)
            c = f.MultEstimators.__getattr__(results_dir_name).V0M.dNpidpT
            sec.add_figure(c)
            c = f.MultEstimators.__getattr__(results_dir_name).V0M.dNKdpT
            sec.add_figure(c)

            # "Backup section with __all__ plots
            sec = latexdoc.add_section(r"Backup")
            sec = latexdoc.add_section(r"$dN/d\eta$")
            for est in settings.considered_ests:
                try:
                    c = f.MultEstimators.__getattr__(results_dir_name).__getattr__(est).__getattr__('dNdeta_summary')
                except AttributeError:
                    # This happens if the estimator does not exist or if
                    # the plot is not present for this estimator
                    continue
                sec.add_figure(c)

            sec = latexdoc.add_section(r"$dN/d\eta$ over MB result")
            for est in settings.considered_ests:
                try:
                    c = f.MultEstimators.__getattr__(results_dir_name).__getattr__(est).__getattr__('dNdeta_MB_ratio_summary')
                except AttributeError:
                    continue
                sec.add_figure(c)

            sec = latexdoc.add_section(r"$P(N_{ch})$ summary of estimators")
            c = f.MultEstimators.__getattr__(results_dir_name).__getattr__('PNch_summary')
            sec.add_figure(c)

            sec = latexdoc.add_section(r"$P(N_{ch}^{\text{estimator 1}})$ binned in $N_{ch}^{\text{estimator 1}}$")
            # The naming convention in this case is a bit confusion. There are essentially two canvas name patterns
            # with the same content. Hence, two canvases
            for est1 in settings.considered_ests:
                for est2 in settings.considered_ests:
                    try:
                        c1 = f.MultEstimators.__getattr__(results_dir_name).__getattr__(est1)\
                            .__getattr__("PNchEst_binned_in_Nch{0}".format(est2))
                        c2 = f.MultEstimators.__getattr__(results_dir_name).__getattr__(est1)\
                            .__getattr__("PNch{0}_binned_in_NchEst".format(est2))
                    except AttributeError:
                        continue
                    for c in [c1, c2]:
                        sec.add_figure(c)

            sec = latexdoc.add_section(r"$\left< p_T \right>$ vs. ref multiplicity")
            for est in settings.considered_ests:
                try:
                    c = f.MultEstimators.__getattr__(results_dir_name).__getattr__(est).__getattr__('mean_pt')
                except AttributeError:
                    continue
                sec.add_figure(c)

            sec = latexdoc.add_section(r"Ratios for various species vs $p_T$")
            for est in settings.considered_ests:
                try:
                    est_dir = f.MultEstimators.__getattr__(results_dir_name).__getattr__(est)
                    canvas_names = [p for p in est_dir.pid_ratios.walk()][-1][-1]
                except AttributeError:
                        continue
                for canvas_name in canvas_names:
                    try:
                        c = est_dir.__getattr__('pid_ratios').__getattr__(canvas_name)
                    except AttributeError:
                        continue
                    sec.add_figure(c)

            sec = latexdoc.add_section(r"Ratios for various species vs ref. multiplicity")
            try:
                ratios_dir = f.MultEstimators.__getattr__(results_dir_name).__getattr__('pid_ratios_vs_refmult')
                canvas_names = [p for p in ratios_dir.walk()][-1][-1]
            except AttributeError:
                    continue
            for canvas_name in canvas_names:
                try:
                    c = ratios_dir.__getattr__(canvas_name)
                except AttributeError:
                    continue
                sec.add_figure(c)

            sec = latexdoc.add_section(r"$dN/dp_T$")
            for est in settings.considered_ests:
                try:
                    c = f.MultEstimators.__getattr__(results_dir_name).__getattr__(est).__getattr__('dNdpT')
                except AttributeError:
                    continue
                sec.add_figure(c)

            sec = latexdoc.add_section(r"$dN_{HM}/dp_T / dN_{MB}/dp_T$")
            for est in settings.considered_ests:
                try:
                    c = f.MultEstimators.__getattr__(results_dir_name)\
                                        .__getattr__(est)\
                                        .__getattr__('pt_hm_div_pt_mb')
                except AttributeError:
                    continue
                sec.add_figure(c)

            sec = latexdoc.add_section(r"$dN_{HM}/dp_T /  dN_{MB}/dp_T \times  \left<N_{MPI}^{MB}\right> / \left<N_{MPI}^{HM}\right>$")
            for est in settings.considered_ests:
                try:
                    c = f.MultEstimators.__getattr__(results_dir_name)\
                                        .__getattr__(est)\
                                        .__getattr__('pt_hm_div_pt_mb_scaled_nMPI')
                except AttributeError:
                    continue
                sec.add_figure(c)

            sec = latexdoc.add_section(r"$nMPI(N_{ch}) summary$")
            c = f.MultEstimators.__getattr__(results_dir_name).__getattr__('nMPI_summary')
            sec.add_figure(c)

            latexdoc.finalize_document()
