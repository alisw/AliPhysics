import re

from rootpy import ROOT
from rootpy.io import root_open

from aligenqa import roofie, utils


def _reduce_single_canvases(results_dirs, gen_names, get_canvas_func, legend_regex):
    fig_comp = roofie.Figure()

    for i, (result_dir, marker) in enumerate(zip(results_dirs, ['opendiamond', 'cross'])):
        c = get_canvas_func(result_dir)
        fig_single = roofie.Figure()
        fig_single.import_plottables_from_canvas(c)
        colors = roofie.get_color_generator('root', ncolors=3)
        for p in fig_single._plottables:
            if legend_regex.match(p['legend_title']):
                color = next(colors)
                fig_comp.add_plottable(p['p'], markerstyle=marker, color=color)
                # set the legend if this is the first file
                if i == 0:
                    fig_comp.add_plottable(None, markerstyle='circle', legend_title=p['legend_title'], color=color)
    # name the lines
    fig_comp.add_plottable(None, markerstyle='opendiamond', legend_title=gen_names[0], color='black')
    fig_comp.add_plottable(None, markerstyle='cross', legend_title=gen_names[1], color='black')
    # configure the plot
    fig_comp.xtitle = fig_single.xtitle
    fig_comp.ytitle = fig_single.ytitle
    fig_comp.legend.title = fig_single.legend.title
    return fig_comp


def compare(args):
    """
    Compare the 'highlight plots' for of two estimators for two given
    triggers. Requires the plots to have been previously prepared by
    running `prepare_plots`.
    """
    ROOT.gROOT.SetBatch(True)
    gen_names = []
    if args.generator_name1:
        gen_names.append(args.generator_name1)
    else:
        gen_names.append(utils.get_generator_name_from_filename(args.input_file1))
    if args.generator_name2:
        gen_names.append(args.generator_name2)
    else:
        gen_names.append(utils.get_generator_name_from_filename(args.input_file2))

    with root_open(args.input_file1, "read") as f1:
        with root_open(args.input_file2, "read") as f2:
            results_dirs = []
            for f, trigger in zip([f1, f2], [args.trigger1, args.trigger2]):
                results_dirs.append(f.MultEstimators.__getattr__("results_post" + trigger))

            latexdoc = roofie.Beamerdoc(author="PWG-MM",
                                        title=r"Comparison of {0} {1} with {2} {3}".format(gen_names[0], args.trigger1,
                                                                                           gen_names[1], args.trigger2),
                                        subtitle="Generator-level QA")

            # Fish the plots we want out of the .root file and if necessary adjust some visual settings
            sec = latexdoc.add_section(r'Comparison of "highlight" plots')

            # filter out an reduced set of plottables
            # to only low mult, high mult and MB.
            # We find those plots by regex-matching the legend
            legend_regex = re.compile(r'.*100'   # contains substring "100"; ie 100% bin edge
                                      r'|^0.0'   # OR starts with "0.0" (highest mult bin edge)
                                      r'|MB')    # OR contains 'MB"

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.EtaLt15.dNdeta_summary, legend_regex)
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.V0M.dNdeta_summary, legend_regex)
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.EtaLt15.dNdeta_MB_ratio_summary, legend_regex)
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.V0M.dNdeta_MB_ratio_summary, legend_regex)
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.EtaLt15.dNdeta_MB_ratio_summary, legend_regex)
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.EtaLt15.PNchEtaLt05_binned_in_NchEst,
                                               legend_regex)
            fig_comp.plot.logy = True
            fig_comp.legend.position = 'tr'
            fig_comp.plot.xmax = 150
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.V0M.PNchEtaLt05_binned_in_NchEst, legend_regex)
            fig_comp.plot.logy = True
            fig_comp.legend.position = 'tr'
            fig_comp.plot.xmax = 150
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.EtaLt15.pid_ratios.proton_over_pich__vs__pt,
                                               legend_regex)
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.EtaLt15.pid_ratios.Xi_over_pich__vs__pt,
                                               legend_regex)
            fig_comp.plot.ymin = 0
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.V0M.pid_ratios.proton_over_pich__vs__pt,
                                               legend_regex)
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.V0M.pid_ratios.Xi_over_pich__vs__pt,
                                               legend_regex)
            fig_comp.plot.ymin = 0
            sec.add_figure(fig_comp)

            # dN/dp_T either for charged particles or individual ones
            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.V0M.dNchdpT,
                                               legend_regex)
            fig_comp.plot.logy = True
            fig_comp.legend.position = 'tr'
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.V0M.dNpdpT,
                                               legend_regex)
            fig_comp.plot.logy = True
            fig_comp.legend.position = 'tr'
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.V0M.dNpidpT,
                                               legend_regex)
            fig_comp.plot.logy = True
            fig_comp.legend.position = 'tr'
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.V0M.dNKdpT,
                                               legend_regex)
            fig_comp.plot.logy = True
            fig_comp.legend.position = 'tr'
            sec.add_figure(fig_comp)

            # only compare a subset of estimators in the following plots; select them via regular expressions
            legend_regex = re.compile(
                r'^\|#eta\|#leq1.5'
                r'|V0M'  # OR V0M
            )
            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: (result_dir.PNch_summary),
                                               legend_regex)
            fig_comp.plot.logy = True
            fig_comp.legend.position = 'tr'
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: (result_dir
                                                                   .pid_ratios_vs_refmult
                                                                   .__getattr__("-2212_2212_div_-211_211")),
                                               legend_regex)
            fig_comp.plot.ymin = 0.04
            fig_comp.plot.ymax = 0.13
            sec.add_figure(fig_comp)

            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: (result_dir
                                                                   .pid_ratios_vs_refmult
                                                                   .__getattr__("3312_div_-211_211")),
                                               legend_regex)
            fig_comp.plot.ymin = 0.0005
            fig_comp.plot.ymax = 0.003
            sec.add_figure(fig_comp)

            # only compare proton, charged kaons and pions
            legend_regex = re.compile(r'^#pi'
                                      r'|K\^{#pm}'  # OR charged kaons
                                      r'|^p$')      # OR proton
            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.EtaLt15.mean_pt, legend_regex)
            sec.add_figure(fig_comp)
            fig_comp = _reduce_single_canvases(results_dirs, gen_names,
                                               lambda result_dir: result_dir.V0M.mean_pt, legend_regex)
            sec.add_figure(fig_comp)

            latexdoc.finalize_document(output_file_name="comparison.tex")
