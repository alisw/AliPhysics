// $Id$

/*! \page README_mchview MUON Tracker visualisation program

A visualisation program, \link mchview.cxx mchview \endlink, is now available to display, in two dimensions
 (3D visu being done within the EVE framework), the tracker chambers. 
 
\section mchview_install Installing the program

\em mchview should be installed together with the rest of AliRoot. One point should be noted though.
\em mchview is using one external file to run properly : a resource file $HOME/.mchviewrc, used
 to "configure" the program and keep some history of the user interaction with it (e.g. the list of recent
 data sources used). If you install a new version of \em mchview, it is recommended to delete this file first, before lauching \em mchview again.
 
\section mchview_navigation Navigating the display

When starting \em mchview (and after the padstore.root file has been created for the first time), 
you'll be presented with two tabs. The first one allow to navigate within the detector, and the 
second one to select data sources to be displayed.
The first tab offers a global view of the 10 tracking chambers. On the right you'll see
 a color palette (used in the display of data sources, see later). On the top (from left to right) are
 navigation buttons (backward and forward), and radio buttons to select the kind of view you'd like : you can
 view things for a given cathode or for a given plane. Note that in some instances those buttons maybe inactive.
On the bottom are three groups of radio buttons, labelled "responder", "outline" and "plot", followed by data source buttons (see later) :
 Each group will contain a number of buttons corresponding to different view levels of the detector (e.g. detection element, manu, buspatch, etc...).
 In \em mchview jargon, what is displayed on screen is called a "painter". The meaning of responder, outline and plot is as follow :
 
 - responder (only one selected at a time) is the type of painter that responds to mouse events. When you mouse over responder painters, they'll
 be highlighted (bold yellow line around them), and some information about them will be displayed in the top right corner.
 If you click on a responder painter, a new view will open to show only this painter. If you meta-click (button2 on a mouse, or alt-click on a Mac for instance)
 on a responder painter, a new view will open, showing the clicked painter and its "dual" (the other cathode or other plane, depending on how the first one
 was defined).
 
 - outline (multiple selection possible) indicates which painter(s) should be outlined. When starting the program, only manus are outlined, but you
 can outline the detection elements, the chambers, etc...
 
  - plot (see later about plotting and data sources) indicates at which level you want to see the data.
   
On the bottom left is a group button used to select which data source should be displayed (empty until you select a data source, see next section).
Next to it will be a list of buttons to select exactly what to plot from a data source (once you've selected a data source only).

Note that whenever you click on a painter and get a new view, you get use the navigation buttons (top left) to go forward and backward, as in 
a web browser for instance. Note also that the \em mchview menu bar contains a "History" menu where you can see (and pick a view) all the views that were opened.

Even before selecting something to plot, at this stage you could use the program to familiarize yourself with the detector structure (aka mapping).

\section mchview_datasource Specifying the data source

The second tab of the \em mchview allows to select one or several data sources. 
Each data source, in turn, will provide one or more "things" to be plotted.
The number of "things" actually depends on the data source.
For instance, a raw data source will allow to plot the mean and the sigma of the pad charges.

Be warned that this part of the program is likely to evolve, as you'll for sure notice that the interface is quite crude for the moment.

From top to bottom, you'll see group of frames used to :

- select from a list of recently used source

- select a raw data source (either by typing in its full pathname, or opening a file dialog). 
The second line in this group is to specify that you want to calibrate the data. Check one of the calibrate buttons, and specify the
location of the OCDB to be used. If that field is not empty (and the corresponding entry is correct, of course), 
the raw data will be calibrated.
The last line in that group is a single check button, to instruct the program to produce histograms of the data (see \ref mchview_histogramming)

- select an OCDB data source (pedestals)

In all the frames, once you've selected or entered the needed information, you'll click on the "Create data source" button,
and a new data source line will appear in the bottom of that tab (and in also in the first tab, that data source will now 
be selectable for plotting). Each data source line indicates the short name of the data source, the full path, and a list of buttons to run, stop, rewind and
remove. Run/Stop/Rewind is only selectable for data sources where the notion of event means something (e.g. for pedestals it won't).
The short name of the data source is as follow :

- RAW# : raw data for run #
- RAW(#) : raw data for simulated run (where run number is always 0, so # here is the number of such data sources opened at the same time)
- HRAW# (or HRAW(#)) : as above, but with histogramming turned on
- (H)CALZ# (or (H)CALZ(#)): as above, but for data where pedestal subtraction has been done (and no gain correction whatsoever)

Note that all the file paths can be local ones or alien ones, if you have a correctly installed alien, and you use a short wrapped to call the \em mchview program.
For instance :
<pre>
alias mchl $HOME/mchview.alien
</pre>

where mchview.alien is a little script :
<pre>
#!/bin/sh
test=`alien-token-info | grep -c expired`
if [ $test -gt 0 ]; then
  echo "Token expired. Getting a new token"
  alien-token-destroy
  alien-token-init
elif [ ! -e /tmp/gclient_env_$UID ]; then
  echo "Getting a token"
  alien-token-init
fi
if [ ! -e /tmp/gclient_env_$UID ]; then
  echo "No token. Exiting"
  exit
fi
source /tmp/gclient_env_$UID
export alien_API_USER=youralienuserid # only needed if different from your local username
mchview $*
</pre>

\section mchview_histogramming Histogramming

Starting at version 0.9 of the \em mchview program, you can now produce histograms of the raw adc values, while running over the
data. For this you have to check the "histogram" button when creating the data source. Please note that turning on the histogram will slow down
a bit the data reading. 
Histograms produced by the program are as compact as possible in order to fit in memory (so they are *not* plain TH1 objects).
Plain TH1 objects are produced later on (on request only), and should be deleted as soon as possible (you have to realize that
1 million TH1 of 4096 channels has no chance to fit in memory...)
Access to the histograms can be done through the GUI, using the right click on any painter.
For extra flexibily, you can also use the root prompt (of the \em mchview program itself).
First get the data object, and then ask the data object to create the histogram(s) you want. Remember to delete those histograms as soon
as you no longer need them :

<pre>
AliMUONPainterRegistry* reg = AliMUONPainterRegistry::Instance();
reg->Print();
AliMUONVTrackerData* data = reg->FindDataSource("HRAW(1)");
TH1* h = data->CreateChannelHisto(707,1025,63);
h->Draw();
delete h;
h = data->CreateManuHisto(707,1025);
etc...
</pre>

You can get histograms for all levels (except PCB) : channel, manu, bus patch, detection element, chamber. See AliMUONVTrackerData doc. for the methods.

\section mchview_savingprinting Saving and printing

From the File menu of the \em mchview application, you can use SaveAs and PrintAs popups to respectively save the current data sources (meaning you can quit
the program and start again with the same filled data sources, without having to rerun on the source) and print the current display.
Printing needs a little bit of polishing (e.g. getting a nice and descriptive title would help a lot), but it's better than nothing.
Note that the \em mchview application now has a \em --use option to reload a previously saved .root file (same effect as using the File/Open menu).

\section mchview_resource Resource file format

The resource file $HOME/.mchviewrc is a normal Root resource file (see TEnv), i.e. a list of "Key:value(s)" lines.

You should avoid to edit it by hand, as most of it is handled by the mchview program itself.

But for information the defined keys so far are :

disableAutoPedCanvas: 1

Use this one to disable the feature that will open automatically 4 canvases each time you open a data source of type "Pedestals".

defaultRange: PED;Mean;0;500|PED;Sigma;0;10|OCC;occ;0;0.01

Use this one to define default ranges for some data sources and their dimensions. In the example above all the pedestals data source will get a display ranging from 0 to 500 for the mean value and 0 to 10 for the sigma value; the occupancy data source (occ dimension) will range from 0 to 0.01.
Those defaults are normally set using, from mchview program (painter master frame tab), using the "Set as default" button, located below the color palette.

In addition, the NumberOfDataSources and DataSource.# (where # ranges from 0 to NumberOfDataSources-1) are used to described the recent opened sources. But those ones should not be edited by hand unless you really know what you are doing.

\section mchview_warnings Important warnings. Please read.

IMPORTANT WARNINGS

In principle, you could have several raw data sources running at the same time. This is NOT currently working. You can have several data sources opened
 at the same time, but not running at the same time. (this has to do with AliRawReader not being thread-safe for the moment).
 
Once you have one or more data sources added, you can go back to first tab and start looking at the data ;-)

This chapter is defined in the READMEmchview.txt file.

*/
