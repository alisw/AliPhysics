import string
import random
import logging
import os

from rootpy import asrootpy, log
from rootpy.plotting import Legend, Canvas, Pad, Graph
from rootpy.plotting.base import Color, MarkerStyle
from rootpy.plotting.utils import get_limits

import ROOT

# from external import husl

# suppress some nonsense logging messages when writing to pdfs.
# Also, setup default logger
log["/ROOT.TCanvas.Print"].setLevel(log.WARNING)
logging.basicConfig(level=logging.DEBUG)
log = log["/roofie"]


def is_plottable(obj):
    """
    Check if the given object is considered a plottable.

    Currently, TH1 and TGraph are considered plottable.
    """
    return isinstance(obj, (ROOT.TH1, ROOT.TGraph))


class Styles(object):
    # Define names of plot layouts:
    class _Default_Style(object):
        pt_per_cm = 28.4527625
        titlefont = 43
        labelfont = 43
        markerSizepx = 4  # number of pixels of the marker

    class Presentation_full(_Default_Style):
        axisTitleSize = 14
        axisLabelSize = 14
        legendSize = 14
        canvasWidth = 340
        canvasHeight = 300
        plot_margins = (.13, .05, .13, .1)   # left, right, bottom, top
        plot_ytitle_offset = 1.15  # factor of the normal offset :P, may lay outside of the canvas

    class Presentation_half(_Default_Style):
        axisTitleSize = 10
        axisLabelSize = 10
        legendSize = 10
        canvasWidth = 170
        canvasHeight = 150
        plot_margins = (.3, .08, .2, .1)
        plot_ytitle_offset = 1

    class Public_full(_Default_Style):
        axisTitleSize = 10
        axisLabelSize = 8
        legendSize = 8
        canvasWidth = 340
        canvasHeight = 300
        plot_margins = (.13, .05, .13, .04)
        plot_ytitle_offset = 1.15


def gen_random_name():
    """Generate a random name for temp hists"""
    return ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(25))


def get_color_generator(palette='root', ncolors=10):
    """
    Returns a generator for n colors.
    Parameters
    ----------
    palette : string
              name of the color palette which should be used
    ncolors : int
              number of colors this palette should have, it might be ignored by some palettes!
    Returns
    -------
    generator :
               colors which can be digested by _rootpy_
    """
    # generated with sns.palplot(sns.color_palette("colorblind", 10))
    if palette == 'colorblind':
        colors = ([(0.0, 0.4470588235294118, 0.6980392156862745),
                   (0.0, 0.6196078431372549, 0.45098039215686275),
                   (0.8352941176470589, 0.3686274509803922, 0.0),
                   (0.8, 0.4745098039215686, 0.6549019607843137),
                   (0.9411764705882353, 0.8941176470588236, 0.25882352941176473),
                   (0.33725490196078434, 0.7058823529411765, 0.9137254901960784)])
    if palette == 'set2':
        colors = ([(0.40000000596046448, 0.7607843279838562, 0.64705884456634521),
                   (0.98131487965583808, 0.55538641635109398, 0.38740485135246722),
                   (0.55432528607985565, 0.62711267120697922, 0.79595541393055635),
                   (0.90311419262605563, 0.54185316071790801, 0.76495195557089413),
                   (0.65371782148585622, 0.84708959004458262, 0.32827375098770734),
                   (0.9986312957370983, 0.85096502233954041, 0.18488274134841617),
                   (0.89573241682613591, 0.76784315109252932, 0.58182240093455595),
                   (0.70196080207824707, 0.70196080207824707, 0.70196080207824707)])
    if palette == 'husl':
        colors = [(0.9677975592919913, 0.44127456009157356, 0.5358103155058701),
                  (0.8616090647292522, 0.536495730113334, 0.19548899031476086),
                  (0.6804189127793346, 0.6151497514677574, 0.19405452111445337),
                  (0.46810256823426105, 0.6699492535792404, 0.1928958739904499),
                  (0.20125317221201128, 0.6907920815379025, 0.47966761189275336),
                  (0.21044753832183283, 0.6773105080456748, 0.6433941168468681),
                  (0.2197995660828324, 0.6625157876850336, 0.7732093159317209),
                  (0.433280341176423, 0.6065273407962815, 0.9585467098271748),
                  (0.8004936186423958, 0.47703363533737203, 0.9579547196007522),
                  (0.962272393509669, 0.3976451968965351, 0.8008274363432775)]
    if palette == 'root':
        # named colors of the ROOT TColor colorwheel are between 800 and 900, +1 to make them look better
        colors = []
        for i in range(0, ncolors):
            colors.append((800 + int(100.0 / ncolors) * i) + 1)
    if colors:
        for color in colors:
            yield color
    else:
        raise ValueError("Unknonw palette")


class Figure(object):
    def __init__(self):
        # User settable parameters:
        self.title = ''
        self.xtitle = ''
        self.ytitle = ''
        self.plot = self.Plot()
        self.legend = self.Legend()

        # Private:
        self._plottables = []
        self.style = Styles.Presentation_full

    class Plot(object):
        logx = False
        logy = False
        gridx = False
        gridy = False
        palette = 'root'
        palette_ncolors = 10
        xmin, xmax, ymin, ymax = None, None, None, None
        frame = None

    class Legend(object):
        title = None
        position = 'tl'

    def _create_legend(self):
        nentries = len([pdic['legend_title'] for pdic in self._plottables if pdic['legend_title'] != ''])
        leg = Legend(nentries, leftmargin=0, rightmargin=0, entrysep=0.01,
                     textsize=self.style.legendSize, textfont=43, margin=0.1, )
        if self.legend.title:
            leg.SetHeader(self.legend.title)
        leg.SetBorderSize(0)  # no box
        leg.SetFillStyle(0)   # transparent background of legend TPave(!)
        return leg

    def _theme_plottable(self, obj):
        try:
            axes = obj.GetXaxis(), obj.GetYaxis()
            for axis in axes:
                axis.SetLabelSize(self.style.axisLabelSize)
                axis.SetLabelFont(self.style.labelfont)
                axis.SetTitleFont(self.style.titlefont)
                axis.SetTitleSize(self.style.axisTitleSize)
            # yaxis only settings:
            axes[1].SetTitleOffset(self.style.plot_ytitle_offset)
        except AttributeError:
            # obj might not be of the right type
            pass
        # apply styles, this might need to get more fine grained
        # markers are avilable in children of TAttMarker
        if isinstance(obj, ROOT.TAttMarker):
            # marker size 1 == 8 px, and never scales with canvas...
            obj.SetMarkerSize(self.style.markerSizepx / 8.0)

    def add_plottable(self, obj, legend_title='', markerstyle='circle', color=None, use_as_frame=None):
        """
        Add a plottable objet to this figure. This function performs a
        copy of the passed object and assigns it a random name. Once
        commited, these should not be touched any more by the user!!!

        Parameters
        ----------
        obj : Hist1D, Graph, None
            A root plottable object; If none, this object will only show up in the legend
        legend_title : string
            Title for this plottable as shown in the legend
        """
        # Make a copy if we got a plottable
        if obj is not None:
            p = asrootpy(obj.Clone(gen_random_name()))
        else:
            p = ROOT.TLegendEntry()
        if isinstance(p, ROOT.TH1):
            p.SetDirectory(0)  # make sure that the hist is not associated with a file anymore!
        self._plottables.append({'p': p,
                                 'legend_title': legend_title,
                                 'markerstyle': markerstyle,
                                 'color': color,
                                 'use_as_frame': use_as_frame,
                                 })

    def import_plottables_from_canvas(self, canvas):
        """
        Import plottables from a canvas which was previously created with roofie

        Parameters
        ----------
        canvas : Canvas
            A canvas which was created with roofie.

        Raises
        ------
        ValueError :
            The given canvas did not have the internal format as expected from roofie canvases
        """
        pad = canvas.FindObject('plot')
        if pad == None:  # "is None" does not work since TObject is not None, but equal to None...
            raise ValueError("Cannot import canvas, since it is not in roofie format.")
        try:
            legend = [p for p in pad.GetListOfPrimitives() if isinstance(p, ROOT.TLegend)][0]
        except IndexError:
            legend_entries = []
        else:
            legend_entries = [e for e in legend.GetListOfPrimitives()]
        # load the plottables but ignore the frame
        plottables = []
        for p in pad.GetListOfPrimitives():
            if is_plottable(p):
                if p.GetName() != "__frame":
                    plottables.append({'p': asrootpy(p.Clone(gen_random_name()))})
                    for legend_entry in legend_entries:
                        if p == legend_entry.GetObject():
                            plottables[-1]['legend_title'] = legend_entry.GetLabel()
                else:
                    self.xtitle = p.GetXaxis().GetTitle()
                    self.ytitle = p.GetYaxis().GetTitle()
        # set legend title if any
        if legend.GetHeader():
            self.legend.title = legend.GetHeader()

        self._plottables += plottables

    def draw_to_canvas(self):
        """
        Draw this figure to a canvas, which is then returned.
        """
        if len(self._plottables) == 0:
            raise IndexError("No plottables defined")
        c = Canvas(width=self.style.canvasWidth,
                   height=self.style.canvasHeight,
                   size_includes_decorations=True)
        if self.legend.position == 'seperate':
            legend_width = .2
            pad_legend = Pad(1 - legend_width, 0, 1., 1., name="legend")
            pad_legend.SetLeftMargin(0.0)
            pad_legend.SetFillStyle(0)  # make this pad transparent
            pad_legend.Draw()
        else:
            legend_width = 0
        pad_plot = Pad(0., 0., 1 - legend_width, 1., name="plot", )
        pad_plot.SetMargin(*self.style.plot_margins)
        pad_plot.Draw()
        pad_plot.cd()

        # awkward hack around a bug in get limits where everything fails if one plottable is shitty...
        xmin, xmax, ymin, ymax = None, None, None, None
        for pdic in self._plottables:
            try:
                limits = get_limits(pdic['p'], logx=self.plot.logx, logy=self.plot.logy)
                # Beware: Python 2 evaluates min/max of None in an undefined way with no error! Wow...
                xmin = min([xmin, limits[0]]) if xmin is not None else limits[0]
                xmax = max([xmax, limits[1]]) if xmax is not None else limits[1]
                ymin = min([ymin, limits[2]]) if ymin is not None else limits[2]
                ymax = max([ymax, limits[3]]) if ymax is not None else limits[3]
            except TypeError:
                # some plottables do not work with this rootpy function (eg. graph without points, tf1)
                # TODO: should be fixed upstream
                pass
        # overwrite these ranges if defaults are given
        if self.plot.xmin is not None:
            xmin = self.plot.xmin
        if self.plot.xmax is not None:
            xmax = self.plot.xmax
        if self.plot.ymax is not None:
            ymax = self.plot.ymax
        if self.plot.ymin is not None:
            ymin = self.plot.ymin

        if not all([val is not None for val in [xmin, xmax, ymin, ymax]]):
            raise TypeError("unable to determine plot axes ranges from the given plottables")

        colors = get_color_generator(self.plot.palette, self.plot.palette_ncolors)

        # draw an empty frame within the given ranges;
        frame_from_plottable = [p for p in self._plottables if p.get('use_as_frame')]
        if len(frame_from_plottable) > 0:
            frame = frame_from_plottable[0]['p'].Clone('__frame')
            frame.Reset()
            frame.SetStats(0)
            frame.xaxis.SetRangeUser(xmin, xmax)
            frame.yaxis.SetRangeUser(ymin, ymax)
            frame.GetXaxis().SetTitle(self.xtitle)
            frame.GetYaxis().SetTitle(self.ytitle)
            self._theme_plottable(frame)
            frame.Draw()
        else:
            frame = Graph()
            frame.SetName("__frame")
            # add a silly point in order to have root draw this frame...
            frame.SetPoint(0, 0, 0)
            frame.GetXaxis().SetLimits(xmin, xmax)
            frame.GetYaxis().SetLimits(ymin, ymax)
            frame.SetMinimum(ymin)
            frame.SetMaximum(ymax)
            frame.GetXaxis().SetTitle(self.xtitle)
            frame.GetYaxis().SetTitle(self.ytitle)
            self._theme_plottable(frame)
            # Draw this frame: 'A' should draw the axis, but does not work if nothing else is drawn.
            # L would draw a line between the points but is seems to do nothing if only one point is present
            # P would also draw that silly point but we don't want that!
            frame.Draw("AL")

        xtick_length = frame.GetXaxis().GetTickLength()
        ytick_length = frame.GetYaxis().GetTickLength()

        for i, pdic in enumerate(self._plottables):
            obj = pdic['p']
            if isinstance(obj, ROOT.TLegendEntry):
                _root_color = Color(pdic['color'])
                _root_markerstyle = MarkerStyle(pdic['markerstyle'])
                obj.SetMarkerStyle(_root_markerstyle('root'))
                obj.SetMarkerColor(_root_color('root'))

            elif isinstance(obj, (ROOT.TH1, ROOT.TGraph, ROOT.TF1)):
                self._theme_plottable(obj)
                obj.SetMarkerStyle(pdic.get('markerstyle', 'circle'))
                if pdic.get('color', None):
                    obj.color = pdic['color']
                else:
                    try:
                        color = next(colors)
                    except StopIteration:
                        log.warning("Ran out of colors; defaulting to black")
                        color = 1
                    obj.color = color
                xaxis = obj.GetXaxis()
                yaxis = obj.GetYaxis()

                # Set the title to the given title:
                obj.title = self.title

                # the xaxis depends on the type of the plottable :P
                if isinstance(obj, ROOT.TGraph):
                    # SetLimit on a TH1 is simply messing up the
                    # lables of the axis to screw over the user, presumably...
                    xaxis.SetLimits(xmin, xmax)
                    yaxis.SetLimits(ymin, ymax)  # for unbinned data
                    # 'P' plots the current marker, 'L' would connect the dots with a simple line
                    # see: https://root.cern.ch/doc/master/classTGraphPainter.html for more draw options
                    drawoption = 'Psame'
                elif isinstance(obj, ROOT.TH1):
                    obj.SetStats(0)
                    xaxis.SetRangeUser(xmin, xmax)
                    yaxis.SetRangeUser(ymin, ymax)
                    drawoption = 'same'
                elif isinstance(obj, ROOT.TF1):
                    # xaxis.SetLimits(xmin, xmax)
                    # yaxis.SetLimits(ymin, ymax)  # for unbinned data
                    drawoption = 'same'
                obj.Draw(drawoption)
            # Its ok if obj is non; then we just add it to the legend.
            else:
                raise TypeError("Un-plottable type given.")
        pad_plot.SetTicks()
        pad_plot.SetLogx(self.plot.logx)
        pad_plot.SetLogy(self.plot.logy)
        pad_plot.SetGridx(self.plot.gridx)
        pad_plot.SetGridy(self.plot.gridy)

        # do we have legend titles?
        if any([pdic.get('legend_title') for pdic in self._plottables]):
            leg = self._create_legend()
            longest_label = 0
            for pdic in self._plottables:
                if not pdic.get('legend_title', False):
                    continue
                leg.AddEntry(pdic['p'], pdic['legend_title'])
                if len(pdic['legend_title']) > longest_label:
                    longest_label = len(pdic['legend_title'])

            # Set the legend position
            # vertical:
            if self.legend.position.startswith('t'):
                leg_hight = leg.y2 - leg.y1
                leg.y2 = 1 - pad_plot.GetTopMargin() - ytick_length
                leg.y1 = leg.y2 - leg_hight
            elif self.legend.position.startswith('b'):
                leg_hight = leg.y2 - leg.y1
                leg.y1 = pad_plot.GetBottomMargin() + ytick_length
                leg.y2 = leg.y1 + leg_hight
            # horizontal:
            if self.legend.position[1:].startswith('l'):
                leg_width = 0.3
                leg.x1 = pad_plot.GetLeftMargin() + xtick_length
                leg.x2 = leg.x1 + leg_width
            elif self.legend.position[1:].startswith('r'):
                leg_width = 0.3
                leg.x2 = 1 - pad_plot.GetRightMargin() - xtick_length
                leg.x1 = leg.x2 - leg_width
            if self.legend.position == 'seperate':
                with pad_legend:
                    leg.Draw()
            else:
                leg.Draw()
        if self.plot.logx:
            pad_plot.SetLogx(True)
        if self.plot.logy:
            pad_plot.SetLogy(True)
        pad_plot.Update()  # needed sometimes with import of canvas. maybe because other "plot" pads exist...
        return c

    def delete_plottables(self):
        """
        Delete all plottables in this figure so that it can be filled with
        new ones while keeping the lables.
        """
        self._plottables = []

    def save_to_root_file(self, in_f, name, path=''):
        """
        Save the current figure to the given root file under the given path
        Parameters
        ----------
        f : TFile
            Root file object open in writable mode
        name : str
            Name for the canvas in the root file
        path : str
            The path where the figure should be saved within the root file
        Returns
        -------
        TFile :
            The file where the object was written to
        """
        f = asrootpy(in_f)
        c = self.draw_to_canvas()
        c.name = name
        try:
            f.mkdir(path, recurse=True)
        except ValueError:
            pass
        f.cd(path)
        success = c.Write()
        if success == 0:
            raise ValueError("Could not write to file!")
        return f

    def save_to_file(self, path, name):
        """
        Save the current figure to the given root file under the given path
        Parameters
        ----------
        name : string
            Name of the file including its extension
        path : string
            Path excluding the file name, relative files are interpreted relative to the working dir
        """
        # check if the name has the right extension
        if len(name.split('.')) != 2:
            raise ValueError("Filename must be given with extension")
        if name.split('.')[1] != 'pdf':
            raise NotImplementedError("Only PDF export is implemented at the moment")
        # strip of tailing / if any
        # this is not compatible with windows, I guess!
        path = path.rstrip('/')
        try:
            os.makedirs(path)
        except OSError:
            pass

        # The order of the following is important! First, set paper
        # size, then draw the canvas and then create the pdf Doin
        # pdf.Range(10, 10) is not sufficient. it just does random
        # sh...
        # Be careful to reset the global gStyle when we are finished. Yeah! Globals!
        # Ok, Root does not like that either...
        # paper_width, paper_height = ROOT.Double(), ROOT.Double()
        # ROOT.gStyle.GetPaperSize(paper_width, paper_height)
        ROOT.gStyle.SetPaperSize(self.style.canvasWidth / self.style.pt_per_cm,
                                 self.style.canvasHeight / self.style.pt_per_cm,)
        c = self.draw_to_canvas()
        c.Print("{0}/{1}".format(path, name))

        # reset the page size
        # ROOT.gStyle.SetPaperSize(paper_width, paper_height)
