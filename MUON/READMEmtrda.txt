// $Id$

/*! 

\page README_mtrda MUON Trigger DA
 
The detector algorithm is implemented for the Muon Trigger in the AliRoot 
framework. The main code is located in MUONTRGda.cxx and it runs in the MUON 
Trigger MON (monitoring).

\section mtrda_s1 The Muon Trigger Calibration

The main goal of the DA is to transfer the configuration files from the detector
data base to the FES and to put them back to the detector data base in order to be
used for the next run, in the case when the analysis of the events in the run finds 
noisy/dead channels different from the ones already notified in the configuration 
files at the start of the run. 
In the current version, the DA will modify only the global crate configuration.

The configuration files stored in the online DB are the following:

- MtgGlobalCrate-[version].dat:   contains the global crate information
- MtgRegionalCrate-[version].dat: contains the regional crate information
- MtgLocalMask-[version].dat:     contains the local mask
- MtgLocalLut-[version].dat:      contains the local LUT 
- MtgCurrent.dat:                 contains the name list of the above files with 
                                  their version and the flag for master/slave
                                  status of the DA
- DAConfig.dat                    configurable parameters for the DA

The copy onto the FES for the modified global masks is done for any value of 
the flag master/slave. The DA creates a file (ExportedFiles.dat) containing the 
name of the files to be transfered by the shuttle. To be able to check the change 
of version of one the files, another file is created containing the last current
list of configuration files: MtgLastCurrent.dat. The Muon trigger electronics can 
run with two types of calibration. The two types of analysis can be done in the same
run containing a mixture of physics events and calibration events.

The SHUTTLE will process the files stored on the FES only in the PHYSICS mode. In STANDALONE
mode the DA will work too, updating the files in the detector data base, but the SHUTTLE will
not process them into the offline data base.

\subsection mtrda_ss1  Dead channels in the global trigger input (with CALIBRATION_EVENT)

The FET pulses are injected to the 21 kchannels.

- The FET mode in MtgGlobalCrate.dat must be set to 0x3
- Data taking (typically 1000 events)
- The DA computes the occupancy of the global input entries, if a channel is not 
responding in more than N% of the events (10% by default), it will be marked as dead 
- The DA updates the global mask file accordingly, adds the file to the data base 
and on the File Exchange Server at the beginning of the next run. 

\subsection mtrda_ss2  Noisy channels in global trigger input (with PHYSICS_EVENT)

This events are used to check the noisy channels triggering at a rate which is 
significantly higher than the expected rate from normal physics events.

- Data taking (typically 1000 events)
- The DA computes the occupancy of the global input entries, if a channel is 
responding in more than N% of the events (10% by default), it will be marked as 
noisy
- The DA updates the global mask file accordingly, adds the file to the data base 
and on the the File Exchange Server at the beginning of the next run. 

\subsection mtrda_ss3 Mixed events (PHYSICS+CALIBRATION)

- If the FET mode in MtgGlobalCrate is not set to 0x3, only PHYSICS events will be used to check for
the noisy channels
- The algorithm is selected according to the event type returned by the raw reader 

\section mtrda_s2 Using the DA Online

With the help of the Control Panel a configuration file is added to the database
(DAConfig.dat) which contains parameters for running the DA (typical values are shown):

- the thresholds for calculating noisy/dead inputs (0.1/0.9)
- the minimum number of events necessary for calculating the input rates (10)
- the maximum number of events to be analyzed in one DA execution (1000000)
- the number of events to skip from the start of run (0)
- the verbosity level of the DA (0, minimum of messages, 1 more messages, 2 print every event)
- enable warnings from the raw data decoder (0, do not show warnings)
- switch between slow/fast payload decoder (0, use the slow decoder)

This file it is not "version"-ed, so it will be not recorded in MtgCurrent.dat.

Test configuration files:

- <a href="http://aliceinfo.cern.ch/static/Offline/dimuon/data/DAConfig.dat">          DAConfig.dat         </a> 
- <a href="http://aliceinfo.cern.ch/static/Offline/dimuon/data/MtgGlobalCrate-1.dat">  MtgGlobalCrate-1.dat </a>
- <a href="http://aliceinfo.cern.ch/static/Offline/dimuon/data/MtgRegionalCrate-1.dat">MtgRegionalCrate-1.dat  </a>
- <a href="http://aliceinfo.cern.ch/static/Offline/dimuon/data/MtgLocalMask-1.dat">    MtgLocalMask-1.dat   </a>
- <a href="http://aliceinfo.cern.ch/static/Offline/dimuon/data/MtgLocalLut-1.dat">     MtgLocalLut-1.dat    </a>
- <a href="http://aliceinfo.cern.ch/static/Offline/dimuon/data/MtgSignature-1.dat">    MtgSignature-1.dat   </a>
- <a href="http://aliceinfo.cern.ch/static/Offline/dimuon/data/MtgCurrent.dat">        MtgCurrent.dat       </a>

\section mtrda_s3 In case of trouble

Please contact:

Franck Manso: manson@clermont.in2p3.fr

or 

Bogdan Vulpescu: vulpescu@clermont.in2p3.fr

This chapter is defined in the READMEmtrda.txt file.
*/

