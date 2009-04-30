/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \file Doxymain.h
/// \brief The main page for code documenation


/*! \mainpage MUON code documentation 

This is the documentation for the MUON simulation and reconstruction code.
It is a mix of general concepts and code implementation details.
It is constantly updated by all dimuon code developers.

The documentation is organized in the thematic pages, defined in the 
README*.txt files, which follow the code organization in the libraries.
Currently there are the documentation pages on
- \ref README_sim
- \ref README_rec
- \ref README_base
- \ref README_mchview 
- \ref README_eve
- \ref README_evaluation 
- \ref README_cosmics
- \ref README_fast
- \ref README_raw 
- \ref README_mapping 
- \ref README_mchda 
- \ref README_mtrda 
- \ref README_calib 
- \ref README_geometry 
- \ref README_trigger 
- \ref README_shuttle 

On this page you will find the first how to run the 
simulation, reconstructin and evaluation chain. More details
and various macros can be found on the other pages.

\section s1 How to check that your aliroot is working well

There is a script file AlirootRun_MUONtest.sh which 
allows for simulating, reconstructing and making the
invariant analysis of the generated Upsilon (1S).
The used configuration file is Config.C in MUON 
directory.

You have to type :
<pre>
$ALICE_ROOT/MUON/AlirootRun_MUONtest.sh [option]
</pre>

The complete list of the option is printed when you call
the script with whatever non valid option, .eg. h:

<pre>
./AlirootRun_MUONtest.sh h
ERROR : extra option not recognized
Usage: AlirootRun_MUONtest.sh options (-SRXsrxn:tg:p:d:c:)
       -S (-s) perform (or not) simulation (default is do it, i.e -S)
       -R (-r) perform (or not) reconstruction (default is do it, i.e. -R)
       -X event (-x) perform (or not) checks and dumps (default is do it for event 5, i.e. -X 5)
       -n nevents (int) number of events to simulate (default 100)
       -t will use OUTDIR as a tmp directory to generate raw data  
       -g seed (uint) seed to be used in simulation (default 1234567)
       -p recoptions (quotified string) reconstruction options to use (default "SAVEDIGITS")
       -d full path to output directory (default /work/projects/alice/trunk/AliRoot/MUON/test_out.100)
       -c full path to configuration file for simulation (default /work/projects/alice/trunk/AliRoot/MUON/Config.C)
</pre>

The results of this test are saved in test_out.nevent directory.
Please note that the CDB (Condition DataBase) is now always *required* 
to perform either simulation or reconstruction. For the moment, a version
of that CDB is stored in CVS, so you should have one already in MUON/Calib
subdirectories.


\section s2 How to check that your aliroot is working VERY well

There is a script file AlirootRun_MUONlongtest.sh which
allows for simulating, reconstructing and making the
-+invariant analysis of the generated Upsilon (1S).
This script generates a large number of Upsilon (20k) 
in order to access differential quantities. 
The used configuration file is Config.C in MUON
directory.

One should really run this script to check if the MUON 
code can process a large number of events WITHOUT errors,
in particular before making important commits !!

You have to type :
<pre>
$ALICE_ROOT/MUON/AlirootRun_MUONtestlong.sh
</pre>
The results of this test are saved in testlong_out/ directory
and will be kept in CVS

(NOTE: the macros performing the calculations/plots MUONefficiency.C 
and MUONplotefficiency.C are also able to handle J/Psi if 
Config.C is modified accordingly )

*/
