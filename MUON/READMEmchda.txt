// $Id$

/*! 

\page README_mchda Tracking DA
 
The detector algorithms are implemented for the Muon Tracking in the AliRoot framework.
We currently have 3 DAs for MCH :

- MUONTRKPEDda.cxx for PEDESTAL runs, running at the end of data taking on each LDC.
- MUONTRKGAINda.cxx for CALIBRATION runs, running at the end of data taking on each LDC.
- MUONTRKOCCda.cxx for PHYSICS runs, running during data taking on each LDC.

\section da_s1 The Muon Tracking Calibration

The Muon tracking chambers needs three types of calibration in order to work properly 
(to be more precise pedestals are required, gains are needed to get the best charge measurement possible, and the occupancy
 is needed in order not to spend all the reconstruction time in hot-spots). 

\subsection da_ss1 Pedestals

The front-end electronics performs an online zero suppression using a threshold level.
Those threshold levels for all channels (~ 1 million) have to be computed in a dedicated
PEDESTALS runs. During this runs the zero suppression is OFF and the pedestal level and the noise is obtained for each channel. The threshold for the FEE is obtained adding the pedestal
level to 3 sigmas of the noise. 

The typical ECS sequence for pedestals is :

- Switch ON the electronics LV
- Boot the CROCUS 
- Configuration 
- Zero suppression OFF
- Data taking (typically 400 events)
- The DA computes the mean and sigma (it runs in each LDC)
- The DA writes one ASCII file per LDC with the results in the File Exchange Server

Then the SHUTTLE process the ASCII files and store the result on the OCDB

\subsection da_ss2 Electronics gain

In order to perform the required spatial resolution or the tracking chambers (~ 100 microns),
we need to calibrate the gain of each channel. The gain is computed using dedicated runs where
a signal is send to the chambers FEE. 

The typical ECS sequence for calibration is :

- Switch ON the electronics LV
- Boot the CROCUS 
- Configuration 
- Zero suppression OFF
- Loop of 10 data taking (typically 400 events) each with a different signal
- The DA computes the mean and sigma (it runs in each LDC) for each run
- The DA write one ASCII file per LDC with the results in the File Exchange Server at the
end of the 10 runs sequence

Then the SHUTTLE process the ASCII files and store the result on the OCDB

\subsection da_ss3 Occupancy

For PHYSICS (or STANDALONE) runs, the MUONTRKOCCda, which is a monitoring DA, keep track of how many times
 each channel has been hit during the run. The output is an ASCII file containing the needed information to 
 compute the occupancy values. This file is written to the DAQ FXS so the SHUTTLE can transfer it to the OCDB.

\section da_s2 Using the DA Online

You have a line command help. To have it just type :

\verbatim
> MUONTRKPEDda.exe -h

******************* ./MUONTRKPEDda.exe usage **********************
./MUONTRKPEDda.exe -options, the available options are :
-h help                   (this screen)

 Input
-f <raw data file>        (default = )

 Output
-a <Flat ASCII file>      (default = )

 Options
-b <output directory>     (default = .)
-c <FES switch>           (default = 1)
-d <print level>          (default = 1)
-g <plot level>           (default = 0)
-i <nb linear points>     (default = 6)
-l <DAC level>            (default = 0)
-m <max date events>      (default = 1000000)
-s <skip events>          (default = 0)
-n <max events>           (default = 1000000)
-r root file data for gain(default = MUONTRKda_gain.data)
-e <execute ped/gain>     (default = ped)
-e <gain create>           make gain & create a new root file
-e <gain>                  make gain & update root file
-e <gain compute>          make gain & compute gains
\endverbatim

 
\section da_s3 Using the DA Offline
 
The DAs normally runs with a RAW data DATE format as input
 If you get a file in root format (e.g. from alien), you can de-rootify it using the
  "deroot" program which is part of aliroot. Note that PED and GAIN DAs work with ROOT input files as well.
  
\section da_s4 In case of trouble

Please contact :

Jean-Luc Charvet : jean-luc.charvet@cea.fr 
or
Alberto Baldisseri : a.baldisseri@cea.fr
or
Laurent Aphecetche : laurent.aphecetche@subatech.in2p3.fr (for OCC DA)

This chapter is defined in the READMEMchda.txt file.
*/

