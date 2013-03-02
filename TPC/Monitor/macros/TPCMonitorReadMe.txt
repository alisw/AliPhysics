////////////// "Usage of the program" ///////////////////////////////////////////////////////////////////////////////

Start of the program:

1. Start aliroot in graphical mode 

2. If monitoring library is not linked to aliroot 
   
       - Load monitoring library from $ALICE_ROOT/lib directory
         gSystem->Load("./lib/tgt_linux/libTPCmon.so")

   Exectute TPCMonitor macro
   .x ${ALICE_ROOT}/Monitor/macros/TPCMonitor.C   	 

// Features of the program /////////////////////////////////////////////////////////////////////////////////////////

1. 'Select Format' Only DATE File and ROOT raw data format are supported by default 
                   If DATE is installed on this machine also online monitoring of raw data is supported
2. 'Select File'   Type filename to be processed.
3. 'Next Event'   (Second last Button): Next EventId to be processed (can be changed to access specific EventId)
                                        The first few events in a file might not be readable --> Go to next event


////////////// "Global view"

Two canvases will pop up showing a global pitcure of the Max adc values for each pad.


Double clicking on one of the sectors in the global canvases will make the corresponding 
IROC and OROC views pop up for the same event.
This can also be achieved by clicking on the corresponding sector button.

'Sector 0' - 'Sector 17'

////////////// "ROC view"

In addition, a pad view for the "Max Adc" for IROC and OROC will pop up.

-Moving with the mouse over a pad in those canvases will make a the 
 time distribution for the pad pop up.

-Double clicking on one of the pads will make a row profile pop up.
 1. Row profile showing Max adc values.
 2. Y-profile for the given position of the pad.
 3  Pad vs. time for all pads in the row.


////////////// "Further options"
 
For a single channel :

'FFT'              : Shows Fourier Transformation.
'Write Channel '   : Writes adc values for channel to ASCII and ROOT (Histogramm) file.

For a sector:

'Show Component'   : Shows only selected components in the "Max Adc" view for a sector.
'Show RMS map'     : Shows pad view for Baseline RMS for IROC and OROC and further distributions.


The Ranges for the calculation of the max adc value and the determination of the baseline can be configured with 

'Conf. Ranges'

The Pedestal calculation is done online by a truncated mean and can be disabled

'Calc BSL (onl)'
'No BSL sub.'

By default a Gamma4 fit is performed to peaks in the time profile for a channel. This can be disabled.

'Disable G4-fit'

For debugging purposes the payload (10bit words format) for each equipment can be written to a file 

'Write 10bit'

The main window size which also determines the canvas sizes can be changed in

${ALICE_ROOT}/TPC/Monitor/AliTPCMonitorConfig.txt
