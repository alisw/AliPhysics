// $Id$

/*! 

\page README_evaluation README evaluation
 

\section evaluation_s1 How to run MUONRecoCheck macro

To check the muon reconstruction by comparing the reconstructed tracks
with the reference tracks made of "AliTrackReference" for the hits in chamber (0..9)
and kinematic informations (TreeK) for the vertex.
This macro can be used to check the track reconstruction e.g. efficiency,
momentum resolution ... but also to make physics analysis whenever
track identification is needed.   

To compile MUONRecoCheck.C
<pre>
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.L $ALICE_ROOT/MUON/MUONRecoCheck.C+
</pre>

To run MUONRecoCheck
<pre>
MUONRecoCheck(nEvent,"geometry.root", "galice.root"); // nEvent = nb of events
</pre>


\section evaluation_s2 Macros for MC studies

For MC studies the classes AliMUONTrackLight and AliMUONPairLight can be 
used in order to fill not only the single muon / dimuon's kinematics (charge, 
pT, rapidity, etc) at the generation AND reconstruction level, but also for 
"decoding" the Pythia output and for the storing of the single muon's history. 
This allows to tag if two muons of a given event come from a certain, well-defined 
process, such as J/psi, Upsilons, correlated open charm or open beauty or the 
low masses or if they are of uncorrelated origin. For open beauty/charm it also 
tags the creation process (pair creation, flavour excitation or gluon splitting). 
The classes also allow to tag feed-down or neutral B meson oscillation and 
has a method that checks whether the reconstructed track is a muon or not.

The macros ReadRecoCocktail.C, DecodeRecoCocktail.C and MergeMuonLight.C 
are examples how to use these two classes. DecodeRecoCocktail.C opens the 
generated files, loops over the events and fills an AliMUONTrackLight object 
for every reconstructed track for which the reference to its generated particle 
could be established, using the AliMUONRecoCheck class. 
It then takes the AliMUONTrackLight objects and forms - event by event - 
AliMUONPairLight objects, on a combinatorial basis. For a given event these 
objects are stored in respective TClonesArrays which are then stored in a tree. 
By default, the produced output file is called "MuonLight.root". 
This root file can then be taken by the macro "ReadRecoCocktail.C" that shows, 
on the example of the reconstructed mass and pT of the AliMUONPairLight object,
how to access the available information. For large statistics, in which many 
individual MuonLight.root files are produced, MergeMuonLight.C can be used 
to merge the files and produce one common output root file.

To read a generation/reconstrution from PDC06 preproduction, and write a file 
with a tree of AliMUONTrackLight / AliMUONPairLight :
go to the directory containing the generation/reconstruction. From there run
aliroot

<pre>
.L DecodeRecoCocktail.C+
DecodeRecoCocktail();
.q
</pre>

To read the file previously generated:
<pre>
aliroot
.L ReadRecoCocktail.C+
ReadRecoCocktail();
.q
</pre>

\section evaluation_s3 Macro MUONCheckDI.C

MUONCheckDI.C performs a consistency check on the methods of the 
AliMUONMCDataInterface and AliMUONDataInterface classes. There are several 
helper methods in these classes which make it easier to fetch data, which 
means there are at least two ways of fetching the data within the same class 
interface. The macro checks to see that the results given by these different 
methods are identical, as they should be.

The macro also inherently exercises the AliMUONMCDataInterface and 
AliMUONDataInterface classes and should be run after any modifications to 
these classes to see if things still work. Putting it another way: 
MUONCheckDI.C is a testing facility for developers of these two classes.

*/
