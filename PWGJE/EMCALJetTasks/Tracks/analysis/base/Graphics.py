#**************************************************************************
#* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
#*                                                                        *
#* Author: The ALICE Off-line Project.                                    *
#* Contributors are mentioned in the code where appropriate.              *
#*                                                                        *
#* Permission to use, copy, modify and distribute this software and its   *
#* documentation strictly for non-commercial purposes is hereby granted   *
#* without fee, provided that the above copyright notice appears in all   *
#* copies and that both the copyright notice and this permission notice   *
#* appear in the supporting documentation. The authors make no claims     *
#* about the suitability of this software for any purpose. It is          *
#* provided "as is" without express or implied warranty.                  *
#**************************************************************************
"""
Graphics module, containing basic ROOT plot helper functionality and
base classes for specific kinds of plots

@author: Markus Fasel
"""

from ROOT import TCanvas,TH1F,TLegend,TPad,TPaveText,TF1
from ROOT import kBlack

class Frame:
    """
    Helper class handling frame drawing in plots
    """
    
    def __init__(self, name, xmin, xmax, ymin, ymax):
        """
        Construct frame with name and ranges for x and y coordinate
        """
        self.__framehist = TH1F(name, "", 100, xmin, xmax)
        self.__framehist.SetStats(False)
        self.__framehist.GetYaxis().SetRangeUser(ymin, ymax)
                        
    def SetXtitle(self, title):
        """
        Set title of the x axis
        """
        self.__framehist.GetXaxis().SetTitle(title)
    
    def SetYtitle(self, title):
        """
        Set title of the y axis
        """
        self.__framehist.GetYaxis().SetTitle(title)
        
    def Draw(self):
        """
        Draw the frame.
        """
        self.__framehist.Draw("axis")

class Style:
    """ 
    Class for plot styles (currently only color and marker)
    """
        
    def __init__(self, color, marker):
        self.__color = color
        self.__marker = marker

    def SetColor(self, color):
        """
        Change color of the graphics object
        """
        self.__color = color

    def SetMarker(self, marker):
        """
        Change marker style of the graphics object
        """
        self.__marker = marker

    def GetColor(self):
        """
        Access color of the graphics object
        """
        return self.__color

    def GetMarker(self):
        """
        Access marker style
        """
        return self.__marker
    
    def DefineROOTPlotObject(self, rootobject):
        """
        Sets the style to the root object
        """
        rootobject.SetMarkerColor(self.__color)
        rootobject.SetMarkerStyle(self.__marker)
        rootobject.SetLineColor(self.__color)
        
class GraphicsObject:
    """
    Container for styled objects
    """
    
    def __init__(self, data, style = None):
        """
        Initialise new graphics object with underlying data (can be TH1 or TGraph(Errors)),
        and optionally a plot style. If no plot style is provided, then the default style (black,
        filled circles) is chosen.
        """
        self.__data = data
        mystyle = Style(kBlack, 20)
        if style:
            mystyle = style
        self.SetStyle(mystyle)
    
    def SetStyle(self, style):
        """
        Initialise underlying object with style
        """
        self.__data.SetLineColor(style.GetColor())
        if not type(self.__data) is TF1:
            self.__data.SetMarkerColor(style.GetColor())
            self.__data.SetMarkerStyle(style.GetMarker())
        
    def GetData(self):
        """
        Provide access to underlying data
        """
        return self.__data
        
    def Draw(self, option = None):
        """
        Draw graphics object. By default, the plot option is 
        "epsame". Option strings will always have the option same
        """
        myoption = "epsame"
        if option:
            myoption = option
            if not "same" in myoption:
                myoption = myoption + "same"
        if type(self.__data) is TF1:
            myoption = "lsame"
        self.__data.Draw(myoption)
    
    def AddToLegend(self, legend, title):
        """
        Add graphics object to a legend provided from outside
        """
        option = "lep"
        if type(self.__data) is TF1:
            option = "l"
        legend.AddEntry(self.__data, title, option)

class PlotBase:
    """
    base class for plot objects
    """
    
    class _FramedPad:
        """
        Defining the pad structure inside the canvas. A pad has a frame with
        axes definition, and optionally a legend and one or several label(s)
        """
        
        class GraphicsEntry:
            """
            Definition of a graphics entry
            """
            
            def __init__(self, graphobject, title = None, addToLegend = False):
                self.__object = graphobject
                self.__title = title
                self.__addToLegend = addToLegend
                
            def __cmp__(self, other):
                """
                Comparison is done accoring to the object title
                """
                # 1st case: either or both of the titles missing
                if not self.__title and not other.GetTitle():
                    return None
                if not self.__title and other.GetTitle():
                    return -1
                if self.__title and not other.GetTitle():
                    return 1
                # second case: both of the titles available
                if self.__title == other.GetTitle():
                    return 0
                if self.__title < other.GetTitle():
                    return -1
                return 1                
                
            def GetObject(self):
                """
                Accessor to graphics object
                """
                return self.__object
            
            def GetTitle(self):
                return self.__title
            
            def IsAddToLegend(self):
                """
                Check whether graphics is foreseen to be added to legend
                """
                return self.__addToLegend
            
            def SetTitle(self, title):
                """
                Change title of the graphics object
                """
                self.__title = title
                
            def SetAddToLegend(self, doAdd):
                self.__addToLegend = doAdd
        
        def __init__(self, pad):
            """
            Constructor, creating a framed pad structure for a TPad
            """
            self.__pad = pad
            self.__Frame = None
            self.__legend = None
            self.__graphicsObjects = []
            self.__labels = []
            
        def DrawFrame(self, frame):
            """
            Draw a frame, defined from outside, within the pad
            The pad becomes owner of the frame
            """
            self.__frame = frame
            self.__frame.Draw()
            
        def DrawGraphicsObject(self, graphics, addToLegend = False, title = None, option = None):
            """
            Draw a graphics object into the pad. If addToLegend is set, then the object is added to to the 
            legend.
            """
            self.__graphicsObjects.append(self.GraphicsEntry(graphics, title, addToLegend))
            graphics.Draw(option)
            
            
        def DefineLegend(self, xmin, ymin, xmax, ymax):
            """
            create a new legend within the frame with the 
            given boundary coordinates
            """
            if not self.__legend:
                self.__legend = TLegend(xmin, ymin, xmax, ymax)
                self.__legend.SetBorderSize(0)
                self.__legend.SetFillStyle(0)
                self.__legend.SetTextFont(42)
                
        def CreateLegend(self, xmin, ymin, xmax, ymax):
            """
            Create Legend from all graphics entries
            """
            if not self.__legend:
                self.DefineLegend(xmin, ymin, xmax, ymax)
                for entry in sorted(self.__graphicsObjects):
                    if entry.IsAddToLegend():
                        self.AddToLegend(entry.GetObject(), entry.GetTitle())
                self.DrawLegend()
            
        def GetLegend(self):
            """
            Provide access to legend
            """
            return self.__legend
        
        def AddToLegend(self, graphicsObject, title):
            """
            Special method adding graphics objects to a legend
            """
            if self.__legend:
                graphicsObject.AddToLegend(self.__legend, title)
            
        def DrawLegend(self):
            """
            Draw the legend
            """
            if self.__legend:
                self.__legend.Draw()
        
        def DrawLabel(self, xmin, ymin, xmax, ymax, text):
            """
            Add a new label to the pad and draw it
            """
            label = TPaveText(xmin, ymin, xmax, ymax, "NDC")
            label.SetBorderSize(0)
            label.SetFillStyle(0)
            label.SetTextFont(42)
            label.AddText(text)
            label.Draw()
            self.__labels.append(label)
        
        def GetPad(self):
            """
            Provide direct access to the pad
            """
            return self.__pad
            
    class _FrameContainer:
        """
        Container for framed pad objects
        """
        
        def __init__(self):
            """
            Create new empty frame container
            """
            self.__Frames = {}
            
        def AddFrame(self, frameID, frame):
            """
            Add a new framed pad to the frame container
            """
            self.__Frames[frameID] = frame
            
        def GetFrame(self, frameID):
            """
            Provide access to frame
            """
            if not self.__Frames.has_key(frameID):
                return None
            return self.__Frames[frameID]

    def __init__(self):
        """
        Initialise new plot
        """
        self._canvas = None
        self._frames = self._FrameContainer()
        
    def _OpenCanvas(self, canvasname, canvastitle, xsize = 1000, ysize = 800):
        """
        Initialise canvas with name, title and sizes
        """
        self._canvas = TCanvas(canvasname, canvastitle, xsize, ysize)
        self._canvas.cd()
        
    def SaveAs(self, filenamebase):
        """
        Save plot to files:
        Creating a file with a common name in the formats
        eps, pdf, jpeg, gif and pdf
        """
        for t in ["eps", "pdf", "jpeg", "gif", "png"]:
            self._canvas.SaveAs("%s.%s" %(filenamebase, t))
            
class SinglePanelPlot(PlotBase):
    
    def __init__(self):
        """
        Initialise single panel plot
        """
        PlotBase.__init__(self)
        
    def _OpenCanvas(self, canvasname, canvastitle):
        """
        Create canvas and add it to the list of framed pads
        """
        PlotBase._OpenCanvas(self, canvasname, canvastitle, 1000, 800)
        self._frames.AddFrame(0, self._FramedPad(self._canvas))
        
    def _GetFramedPad(self):
        """
        Access to framed pad
        """
        return self._frames.GetFrame(0)
    
class MultipanelPlot(PlotBase):
    """
    Base Class For multiple panel plots
    """
    
    def __init__(self, nrow, ncol):
        """
        Create new Multi-panel plot with a given number of row and cols
        """
        PlotBase.__init__(self)
        self.__nrow = nrow
        self.__ncol = ncol
        
    def _OpenCanvas(self, canvasname, canvastitle, xsize, ysize):
        """
        Create new canvas and split it into the amount of pads as defined
        """
        PlotBase._OpenCanvas(self, canvasname, canvastitle, xsize, ysize)
        self._canvas.Divide(self.__ncol, self.__nrow)
        
    def _OpenPad(self, padID):
        """
        Create new framed pad in a multi-panel plot for a given pad ID
        """
        if padID < 0 or padID > self.__GetMaxPadID():
            return None
        mypad = self._GetPad(padID)
        if not mypad:
            mypad = self._FramedPad(self._canvas.cd(padID+1))
            self._frames.AddFrame(padID, mypad)
        return mypad
    
    def _OpenPadByRowCol(self, row, col):
        """
        Create new framed pad in a multi-panel plot for a given row an col
        """
        return self._OpenPad(self.__GetPadID(row, col))
    
    def _GetPad(self, padID):
        """
        Access to Pads by pad ID
        """
        return self._frames.GetFrame(padID)
    
    def _GetPadByRowCol(self, row, col):
        """
        Access Pad by row and col
        """
        return self._frames.GetFrame(self.__GetPadID(row, col))
    
    def __GetPadID(self, row, col):
        """
        Calculate ID of the pad
        """
        if (row < 0 or row >= self.__nrow) or (col < 0 or col >= self.__ncol):
            return -1
        return 1 + row * self.__ncol + col
    
    def __GetMaxPadID(self):
        """
        Calculate the maximum allowed pad ID
        """
        return 1 + self.__ncol * self.__nrow
    
class TwoPanelPlot(MultipanelPlot):
    """
    A plot with two panels
    """
    
    def __init__(self):
        """
        Initialise two-panel plot
        """
        MultipanelPlot.__init__(self, 1, 2)
        
    def _CreateCanvas(self, canvasname, canvastitle):
        """
        Create Canvas with the dimensions of a four-panel plot
        """
        MultipanelPlot._OpenCanvas(self, canvasname, canvastitle, 1000, 500)
    
class FourPanelPlot(MultipanelPlot):
    """
    A plot with four (2x2) panels
    """
    
    def __init__(self):
        """
        Initialise four-panel plot
        """
        MultipanelPlot.__init__(self, 2, 2)
        
    def _OpenCanvas(self, canvasname, canvastitle):
        """
        Create Canvas with the dimensions of a four-panel plot
        """
        MultipanelPlot._OpenCanvas(self, canvasname, canvastitle, 1000, 1000)