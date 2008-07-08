// $Id$

/*! 

\page README_mtrda Trigger DA
 
The detector algorithm is implemented for the Muon Trigger in the AliRoot framework.
The main code is located in MUONTRGda.cxx and it runs in the MUON Trigger LDC.

\section da_s1 The Muon Trigger Calibration

The main goal of the DA is the transfert of the modified configuration files to the FES.
The configuration files stored in the online DB are the following:

- MtgGlobalCrate-<version>.dat:   contains the global crate information
- MtgRegionalCrate-<version>.dat: contains the regional crate information
- MtgLocalMask-<version>.dat:     contains the local mask
- MtgLocalLut-<version>.dat:      contains the local LUT 
- MtgCurrent.dat:                 contains the name list of the above files with their version 
                                and the flag for master/slave status on the DA

The copy onto the FES for the modified local masks is only done when the flag is set to master for the DA.
The DA creates a file (ExportedFiles.dat) containing the name of the files to be transfert by the shuttle.
To be able to check the change of version of one the files, another file is created containing the last current
list of configuration files: MtgLastCurrent.dat.
The Muon trigger electronics could run with two types of calibration:

\subsection da_ss1  ELECTRONICS_CALIBRATION_RUN (calib)

This procedure allows to check dead channels. The FET pulses are injected to the 21 kchannels.

The typical ECS sequence for calib is :

- Switch ON the electronics LV
- Load Configuration via the MTS package
- Enable FET pulse 
- Data taking (typically 100 events)
- The DA computes the occupancy, if a channel is not responding in 80% of the case, it will be marked as dead 
- The DA update the local mask file accordingly and put it on the File Exchange Server

Then the SHUTTLE process the ASCII files and store the configuration on the OCDB

This option is disable for the moment.

\subsection da_ss2  DETECTOR_CALIBRATION_RUN (ped)

This procedure checks the noisy channels. A normal physics run is performed.

The typical ECS sequence for calibration is :

- Switch ON the electronics LV
- Load Configuration via the MTS package
- Data taking (typically 100 events)
- The DA computes the occupancy, if a channel is responding in 80% of the case, it will be marked as noisy
- The DA update the local mask file accordingly and put it on the File Exchange Server

Then the SHUTTLE process the ASCII files and store the configuration on the OCDB

\section da_s2 Using the DA Online

You have a line command help. To have it just type :

\verbatim
> MUONTRGda.exe -h

******************* MUONTRGda.exe usage **********************
MUONTRGda.exe -options, the available options are :
-h help                   (this screen)

 Input
-f <raw data file>        (default = )

 output
-r <root file>            (default = none)

 Options
-t <threshold values>     (default = 0.2)
-d <print level>          (default = 0)
-s <skip events>          (default = 0)
-n <max events>           (default = 1000000)
-e <execute ped/calib>    (default = nil)

\endverbatim

 
\section da_s3 In case of trouble

Please contact:

Franck Manso: manson@clermont.in2p3.fr

or 

Christian Finck: finck@subatech.in2p3.fr (until 01 Septembre 08)


This chapter is defined in the READMEmtrda.txt file.
*/

