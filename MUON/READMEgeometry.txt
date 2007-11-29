// $Id$

/*! 

\page README_geometry README geometry


\section geometry_s1 General Information about MUON Geometry

Our geometry is described in the geometry builder classes.
Main geometrical constants are set in the class AliMUONConstants.
The code can then generate the geometry data files
transform.dat and svmap.dat (see description below) via the macro  
MUONGenerateGeometryData.C (more info below).

The geometry data files have to be recreated each time the code 
of the geometry is modified. The info (well updated) in this files 
(svmap) is need during the simulation.
We can also decide to use the transform.dat file as input of our 
geometry. This allows for changing the position of our detection elements
and/or half-planes (half-chambers in code jargon) without modifying 
and recompiling the code. 

First step in the official aliroot simulation process is to create 
the geometry.root file from the builders to build the MUON geometry 
within the geometrical modeler framework of root. 
Then  aliroot takes the geometry.root file as a unique geometrical 
info of our apparatus during the generation and the reconstruction
and analysis (if needed)

Misalignments are in the official AliRoot code applied to the geometry.root
file.


\section geometry_s2 How to check the Geometry with the new Geometrical modeler

\see ftp://root.cern.ch/root/doc/chapter16.pdf
\see http://agenda.cern.ch/fullAgenda.php?ida=a05212

<pre>
gAlice->Init("$ALICE_ROOT/MUON/Config.C");
gGeoManager->GetMasterVolume()->Draw();
</pre>


\section geometry_s3  How to check the overlap with the new Geometrical modeler

\see  ftp://root.cern.ch/root/doc/chapter16.pdf
\see  http://agenda.cern.ch/fullAgenda.php?ida=a05212

<pre>
gAlice->Init("$ALICE_ROOT/MUON/Config.C");
gGeoManager->CheckOverlaps();
gGeoManager->PrintOverlaps();
</pre>


\section geometry_s4 Macro  MUONGenerateGeometryData.C
						
Macro for generating the geometry data files

Geometry data files:
- MUON/data/transform.dat file contains the transformations
data (translation and rotation) for all alignable objects
(modules & detection elements)
- MUON/data/svmap.dat file contains all the information to link 
each geant volume (it can be extended to other virtual MC) with
a detection element. The point here is that a given detection
element, i.e. a slat chamber can consist of more geant volumes.
the correspondence is then defined in an input file.
Each time there is a change in the definition of MC geometry, these
input files must be re-generated via the macro  
MUONGenerateGeometryData.C

To be run from aliroot:
<pre>
.x MUONGenerateGeometryData.C
</pre>

The generated files do not replace the existing ones
but have different names (with extension ".out").
Replacement with new files has to be done manually.


\section geometry_s5 Macros to generate Mis-alignment data
						
Macros for generating the geometry mis-alignment data: 
- MakeMUONFullMisAlignment.C
- MakeMUONResMisAlignment.C
- MakeMUONZeroMisAlignment.C

To be run from aliroot:
<pre>
.x MakeMUONFullMisAlignment.C etc.
</pre>

If the environment variable TOCDB is not set to "kTRUE",
the misalignment data are generated in a local file:
(MUONFullMisalignment.root, etc.)

If the data are stored in CDB, the storage can be specified in 
the environment variable STORAGE. The misalignment data are then
generated in the CDB folder (defaults are ResMisAlignCDB and FullMisAlignCDB
in the working directory). Inside the local CDB the path for the
alignment data is (and must be) "MUON/Align/Data/".
Residual misalignment: Default is our current estimate of
misalignment after all our alignment procedure has been applied.
Full misalignment: Default is our current estimate of initial
misalignment.


\section geometry_s6 How to check the alignment software

The script AlirootRun_MUONtestAlign.sh  allows you to check the software for
the alignment with physics tracks. The script will:
- Generate a misaligned geometry in a local CDB (default FullMisAlignCDB)
- Simulate 1000 events using previously misaligned geometry
- Reconstruct the events using perfect geometry
- Run the alignment code over the above events using MUONAlignment.C

To run you need to type:
<pre>
$ALICE_ROOT/MUON/AlirootRun_MUONtestAlign.sh
</pre>

The results of the test are saved in test_align/ directory. The file measShifts.root
contains useful graphs for studying the alignment performances. A local CDB
containing the realigned geometry is also created (default is ReAlignCDB). The
file $ALICE_ROOT/MUON/data/transform2ReAlign.dat contains the
transformations describing the realigned geometry to be compared with the
used misaligned geometry $ALICE_ROOT/MUON/data/transform2.dat.

IMPORTANT NOTE: For a useful test of the alignment performances, the
order of 100 000 tracks is needed, it is then advisable to generate and
reconstruct enough events separately and run MUONAlignment.C providing a file list
afterwards.

\section geometry_s7 Macro MUONCheckMisAligner.C

The macro MUONCheckMisAligner.C performs the misalignment on an existing muon 
arm geometry based on the standard definition of the detector elements.

It uses AliMUONGeometryAligner : 
- creates a new AliMUONGeometryTransformer and AliMUONGeometryAligner
- reads the transformations in from the transform.dat file (make sure that
this file is the _standard_ one by comparing it to the one in CVS)
- creates a second AliMUONGeometryTransformer by misaligning the existing 
one using AliMUONAligner::MisAlign

User has to specify the magnitude of the alignments, in the Cartesian 
co-ordiantes (which are used to apply translation misalignments) and in the
spherical co-ordinates (which are used to apply angular displacements)

User can also set misalignment ranges by hand using the methods : 
SetMaxCartMisAlig, SetMaxAngMisAlig, SetXYAngMisAligFactor
(last method takes account of the fact that the misalingment is greatest in 
the XY plane, since the detection elements are fixed to a support structure
in this plane. Misalignments in the XZ and YZ plane will be very small 
compared to those in the XY plane, which are small already - of the order 
of microns)

Default behavior generates a "residual" misalignment using gaussian
distributions. Uniform distributions can still be used, see 
AliMUONGeometryAligner.

User can also generate module misalignments using SetModuleCartMisAlig
and SetModuleAngMisAlig
Note : If the detection elements are allowed to be misaligned in all
directions, this has consequences for the alignment algorithm, which 
needs to know the number of free parameters. Eric only allowed 3 : 
x,y,theta_xy, but in principle z and the other two angles are alignable
as well.  


\section geometry_s8 Geometry data files description
 
\subsection geometry_s8_sub1 transform.dat
 
 List of transformations for chambers geometry modules and detection
 elements; in format:
<pre> 
 KEY   ID  [nofDE]  pos: posX posY posZ  rot: theX phiX theY phiY theZ phiZ
  
 where  KEY  = CH or DE
        ID   = chamberId or detElemId
	pos: posX posY posZ  = position in cm
	rot: theX phiX theY phiY theZ phiZ = rotation angles as in Geant3 in deg
</pre>

\subsection geometry_s8_sub2  svmap.dat

 Map of sensitive volumes to detction element Ids;
 in format:

<pre> 
 KEY  volpath  detElemId
  
 where  KEY  = SV
	volpath   = volume path in format /volname1_copyNo1/volname2_copyNo2/...
        detElemId = detection element Id
 </pre>


*/
