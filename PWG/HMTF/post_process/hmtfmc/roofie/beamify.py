"""
A class mirroring a subset of the Latex beamer functinality.

It is intended as an easy way to compile many roofie figures with small annotations into a valid tex document.
At the same time, it tries to abstract earthly problems like uniqueness of the filenames of the plots.
This is at the expense of customizability of the folder structure
"""

import os
import random
import re
import string
import subprocess
import textwrap

from .figure import Figure

import ROOT.TCanvas

try:
    from subprocess import check_output
except ImportError:
    # check_output not avialable in python 2.6
    def check_output(*popenargs, **kwargs):
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd)
        return output


def _safe_folder_name(unsafe):
    t = re.compile("[a-zA-Z0-9.,_-]")
    return "".join([ch for ch in unsafe if t.match(ch)])


class Beamerdoc(object):
    def __init__(self, author, title):
        self.title = title
        self.author = author
        self.sections = []
        self.output_dir = './{0}/'.format(_safe_folder_name(self.title))
        self.preamble = textwrap.dedent(r"""
        \documentclass[xcolor=dvipsnames]{{beamer}}

        \usepackage{{graphicx,subfigure,url, tikz}}

        % example themes
        % \usetheme[nat,style=simple]{{Frederiksberg}}

        % put page numbers into the footer
        \setbeamertemplate{{footline}}[frame number]{{}}
        % remove navigation symbols
        \setbeamertemplate{{navigation symbols}}{{}}
        \author{{{author}}}
        \institute[NBI, Copenhagen]{{Niels Bohr Institute, Copenhagen\\HMTF MC benchmark studies\\Supervisor: Michele Floris}}

        \title{{{title}}}
        \begin{{document}}

        \frame[plain]{{\titlepage}}
        \frame{{\tableofcontents}}

        """)
        self.postamble = r"""\end{document}\n"""

    class Section(object):
        def __init__(self, document, sec_title):
            self.document = document
            self.title = sec_title
            self.figures = []
            self.frames = []
            self.frame_template = textwrap.dedent(r"""
            \begin{{frame}}
            \frametitle{{{title}}}
            \begin{{columns}}
            \begin{{column}}{{.45\textwidth}}
            {0}\\
            {1}
            \end{{column}}
            \begin{{column}}{{.45\textwidth}}
            {2}\\
            {3}
            \end{{column}}
            \end{{columns}}
            \end{{frame}}

            """)

        def add_figure(self, fig):
            self.figures.append(fig)

        def _make_section_body(self):
            """
            Convert this section to latex. This writes the figures to disc

            Returns
            -------
            string :
                Latex code for this section with linkes to the figures already included
            """
            fig_paths = self._write_figures_to_disc()
            section_body = '\section{{{0}}}'.format(self.title)
            figs_per_frame = 4
            for frame_num in xrange(0, len(fig_paths), figs_per_frame):
                ig_cmds = ['' for i in range(figs_per_frame)]
                for fig_num_frame, path in enumerate(fig_paths[frame_num:(frame_num + figs_per_frame)]):
                    ig_cmds[fig_num_frame] = format(r'\includegraphics[width=\textwidth]{{{0}}}'.format(path))
                section_body += self.frame_template.format(title=self.title, *ig_cmds)
            return section_body

        def _write_figures_to_disc(self):
            """
            Write the figures of this section to disc.

            Returns
            -------
            list :
                Paths to the files written to disc, relative to where the tex file will be
            """
            paths = []
            fig_folder_safe = os.path.join('figures/', _safe_folder_name(self.title))
            for fig in self.figures:
                rand_name = ''.join(random.choice(string.ascii_letters) for _ in range(5)) + ".pdf"
                fig_path_from_latex_root = os.path.join(self.document.output_dir, fig_folder_safe)
                if isinstance(fig, Figure):
                    fig.save_to_file(fig_path_from_latex_root, rand_name)
                if isinstance(fig, ROOT.TCanvas):
                    # make sure the folder exists
                    try:
                        os.makedirs(fig_path_from_latex_root)
                    except OSError:
                        pass
                    # root needs to draw the canvas first otherwise it crashes, sometimes...
                    fig.Draw()
                    fig.SaveAs(os.path.join(fig_path_from_latex_root, rand_name))
                paths.append(os.path.join("./", fig_folder_safe, rand_name))
            return paths

    def add_section(self, sec_title):
        self.sections.append(self.Section(self, sec_title=sec_title))
        return self.sections[-1]

    def _create_latex_and_save_figures(self):
        body = ''
        body += self.preamble.format(title=self.title, author=self.author)
        for sec in self.sections:
            body += sec._make_section_body()
        body += self.postamble
        return body

    def finalize_document(self, output_file_name="summary.tex"):
        """
        Assamble the latex document and run the compile command
        """
        try:
            os.makedirs(self.output_dir)
        except OSError:
            pass

        with file(os.path.abspath(self.output_dir) + "/" + output_file_name, 'w')as f:
            f.write(self._create_latex_and_save_figures())

        cmd = ['pdflatex', '-file-line-error', '-interaction=nonstopmode', format(output_file_name)]
        check_output(cmd, cwd=os.path.abspath(self.output_dir))
        try:
            # only cache errors for the second compiling try, since the first  might complain about old stuff
            check_output(cmd, cwd=os.path.abspath(self.output_dir))
        except subprocess.CalledProcessError:
            print "An error occured while compiling the latex document. See 'summary.log' for details"

# from roofie.figure import Figure
# from rootpy.plotting import Hist1D
# f = Figure()
# h = Hist1D(10, 0, 10)
# h.fill(5)
# f.add_plottable(h)
# latexdoc = Beamerdoc()
# latexdoc.author = "Christian Bourjau"
# latexdoc.title = "Test"
# sec = latexdoc.add_section("My section")
# sec.add_figure(f)
# latexdoc.finalize_document()
