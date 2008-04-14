// $Id$

/*! 

\page README_mchda Tracking DA
 
The detector algorithm is implemented for the Muon Tracking in the AliRoot framework.
The main code is located in MUONTRKda.cxx and it runs in each Muon Tracking LDC.
 
\section da_s1 The Muon Tracking Calibration

The Muon tracking chambers needs two types of calibration in order to work properly :

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

\section da_s2 Using the DA Online

You have a line command help. To have it just type :

\verbatim
> MUONTRKda.exe -h

******************* ./MUONTRKda.exe usage **********************
./MUONTRKda.exe -options, the available options are :
-h help                   (this screen)

 Input
-f <raw data file>        (default = )

 Output
-a <Flat ASCII file>      (default = )

 Options
-b <output directory>     (default = .)
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
 
For the time being the DA can be used only with a RAW data DATE format as input.
The development of an offline version is under way.

\section da_s4 In case of trouble

Please contact :

Jean-Luc Charvet : jean-luc.charvet@cea.fr 
or
Alberto Baldisseri : a.baldisseri@cea.fr


This chapter is defined in the READMEMchda.txt file.
*/

