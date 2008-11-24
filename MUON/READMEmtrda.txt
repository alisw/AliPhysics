// $Id$

/*! 

\page README_mtrda Trigger DA
 
The detector algorithm is implemented for the Muon Trigger in the AliRoot 
framework. The main code is located in MUONTRGda.cxx and it runs in the MUON 
Trigger MON (monitoring).

\section da_s1 The Muon Trigger Calibration

The main goal of the DA is the transfer of the modified configuration files to 
the FES and to put them in the detector data base. In the current version, the 
DA will modify only the global crate configuration.

The configuration files stored in the online DB are the following:

- MtgGlobalCrate-<version>.dat:   contains the global crate information
- MtgRegionalCrate-<version>.dat: contains the regional crate information
- MtgLocalMask-<version>.dat:     contains the local mask
- MtgLocalLut-<version>.dat:      contains the local LUT 
- MtgCurrent.dat:                 contains the name list of the above files with their version 
                                  and the flag for master/slave status on the DA

The copy onto the FES for the modified global masks is done for any value of 
the flag master/slave. The DA creates a file (ExportedFiles.dat) containing the 
name of the files to be transfered by the shuttle. To be able to check the change 
of version of one the files, another file is created containing the last current
list of configuration files: MtgLastCurrent.dat. The Muon trigger electronics can 
run with two types of calibration. New: the two types can be done in the same
run containing a mixture of physics events with calibration events injected every 
50 seconds.

\subsection da_ss1  ELECTRONICS_CALIBRATION_RUN (calibration)

This procedure allows to check dead channels. The FET pulses are injected to the 21 kchannels.

The typical ECS sequence for calib is :

- Switch ON the electronics LV
- Load Configuration via the MTS package
- Enable FET pulse 
- Data taking (typically 1000 events)
- The DA computes the occupancy of the global input entries, if a channel is not 
responding in more than N% of the events (10% by default), it will be marked as dead 
- The DA updates the global mask file accordingly, adds the file to the data base 
and on the File Exchange Server at the beginning of the next run. 

Then the SHUTTLE process the ASCII files and store the configuration on the OCDB.

\subsection da_ss2  DETECTOR_CALIBRATION_RUN (pedestal)

This procedure checks the noisy channels. Normal physics events are used.

The typical ECS sequence for calibration is :

- Switch ON the electronics LV
- Load Configuration via the MTS package
- Data taking (typically 1000 events)
- The DA computes the occupancy of the global input entries, if a channel is 
responding in more than N% of the events (10% by default), it will be marked as 
noisy
- The DA updates the global mask file accordingly, adds the file to the data base 
and on the the File Exchange Server at the beginning of the next run. 

Then the SHUTTLE process the ASCII files and store the configuration on the OCDB.

\section da_s2 Using the DA Online

With the help of the Control Panel a configuration file is added to the database
(DAConfig.txt) which contains parameters for running the DA:

- the thresholds for calculating noisy/dead inputs
- the minimum number of events necessary for calculating the input rates
- the maximum number of events to be analyzed in one DA execution
- the number of events to skip from the start of run
- the verbosity level of the DA
- enable warnings from the raw data decoder

This file it is not "version"-ed, so it will be not recorded in MtgCurrent.dat.
 
\section da_s3 In case of trouble

Please contact:

Franck Manso: manson@clermont.in2p3.fr

or 

Bogdan Vulpescu: vulpescu@clermont.in2p3.fr

This chapter is defined in the READMEmtrda.txt file.
*/

