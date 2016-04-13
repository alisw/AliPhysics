/*! \page READMEemcfw The EMCAL Analysis framework 

The EMCAL framework consists of three parts:

- \subpage READMEcorefw
- \subpage READMEtrgfw
- \subpage READMEfwtasks

Goal of the EMCAL framework is to provide the user with all functionalities
to steer the different steps in the preparation chain, handle data, process
trigger quantities, and combine common functionalities. Furthermore helper
classes for various issues are provided.

\section s1 How to run it

The analysis chain has to be configured with several framework wagons. They
are combined by several Add macros.

\subsection c1 The setup task wagon

The setup task wagon combines handles the connection to OCDB, OADB and initializes
the geometry, used in further wagons. One can add it to the train by

<pre>
gSystem->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
AddTaskEmcalSetup();
</pre>

Paths to OCDB and OADB can be handled by the macro paramters.

\subsection c2 The EMCAL preparation wagon

Recalibration and clusterization is handled by the EMCAL preparation wagon. 

<pre>
gSystem->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPreparation.C");
AddTaskEmcalPreparation();
</pre>

*For suggestions and fixes, write to ...*
*/
