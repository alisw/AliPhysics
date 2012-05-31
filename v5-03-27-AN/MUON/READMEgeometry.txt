// $Id$

/*! 

\page README_geometry Geometry


\section geometry_s1 General Information about MUON Geometry

Our geometry is described in the geometry builder classes.
The main geometrical constants are set in the class AliMUONConstants.
The geometry is built from the code during running simulation
and it is automatically exported in a geometry.root file
via the framework. Then  aliroot takes this geometry.root file as 
a unique geometrical info of our apparatus during the generation 
and the reconstruction and analysis (if needed)

The macros MakeMUONZeroMisAlignment.C, MakeMUONResMisAlignment.C
and MakeMUONFullMisAlignment.C generate the mis-alignment
data (see more in the chapter \ref geometry_s4 below).

The code can also generate the special geometry 
data files, transform.dat and svmap.dat, via the macro  
MUONGenerateGeometryData.C (see more in the chapter \ref geometry_s5 below).
The svmap.dat data file have to be recreated each time the code 
of the geometry is modified. The info (well updated) in this file 
is needed during the simulation.
We can also decide to use the transform.dat file as input of our 
geometry. This allows for changing the position of our detection elements
and/or half-planes (half-chambers in code jargon) without modifying 
and recompiling the code. 

Misalignments are in the official AliRoot code applied to the geometry.root
file.


\section geometry_s2 How to check the geometry with the Root geometrical modeler

\see ftp://root.cern.ch/root/doc/chapter16.pdf
\see http://agenda.cern.ch/fullAgenda.php?ida=a05212

<pre>
TGeoManager::Import("geometry.root");
gGeoManager->GetMasterVolume()->Draw();
</pre>

A helper macro for adding and removing volumes in the
scene, MUONGeometryViewingHelper.C is also available.
 

\section geometry_s3  How to check the overlaps with the Root geometrical modeler

\see  ftp://root.cern.ch/root/doc/chapter16.pdf
\see  http://agenda.cern.ch/fullAgenda.php?ida=a05212

<pre>
TGeoManager::Import("geometry.root");
gGeoManager->CheckOverlaps(0.001);
gGeoManager->PrintOverlaps();
</pre>

More extensive, but also more time consuming checking,
can be performed in this way:
<pre>
gGeoManager->CheckGeometryFull(1000000,0,0,0,"o"); >& check_full.out
</pre>
Then, you will find in the output file \em check_full.out the list of
volumes where any overlaps have been detected. As TGeoManager
does not remember all overlaps found during checking,
in order to investigate them, one has to re-run the checking for 
each listed volume:
<pre>
gGeoManager->FindVolumeFast("MyVolume")->CheckOverlaps(0.001, "s");
gGeoManager->PrintOverlaps(); >& overlaps_MyVolume.txt 
</pre>
At this stage the overlaps found for the selected volume can be also browsed 
with TBrowser. Sometimes it happens that the reported overlapping
volumes are assemblies and nothing is visualized on the scene
when clicking on the overlap icon in the browser.
In this case you can use the function setDaughtersVisibility()
from the MUONGeometryViewingHelper.C macro, which propagates the
visibility setting through all assembly levels up to the real
volumes.

\section geometry_s4 Macro  MUONGenerateGeometryData.C
						
Macro for generating the geometry data files:
- MUON/data/svmap.dat file contains all the information to link 
each geant volume (it can be extended to other virtual MC) with
a detection element. The point here is that a given detection
element, i.e. a slat chamber can consist of more geant volumes.
the correspondence is then defined in an input file.
Each time there is a change in the definition of MC geometry, these
input files must be re-generated via the macro  
MUONGenerateGeometryData.C
- MUON/data/transform.dat file contains the transformations
data (translation and rotation) for all alignable objects
(modules & detection elements)

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
.x MakeMUONFullMisAlignment.C
</pre>

etc.

If the environment variable TOCDB is not set to "kTRUE",
the misalignment data are generated in a local file:
MUONfullMisalignment.root, etc.

If the data are stored in CDB, the storage can be specified in 
the environment variable STORAGE. The misalignment data are then
generated in the CDB folder (defaults are ResMisAlignCDB and FullMisAlignCDB
in the working directory). Inside the local CDB the path for the
alignment data is (and must be) "MUON/Align/Data/".
Residual misalignment: Default is our current estimate of
misalignment after all our alignment procedure has been applied.
Full misalignment: Default is our current estimate of initial
misalignment.

The mis-alignment data can be then retrieved from a file
and applied to ideal geometry in this way.

<pre>
TGeoManager::Import("geometry.root");
TFile f("MUONfullMisalignment.root"); 
TClonesArray* misAlignObjsArray = (TClonesArray*)f.Get("MUONAlignObjs");
AliGeomManager::ApplyAlignObjsToGeom(*misAlignObjsArray);
</pre>

Mis-aligned geometry can be then inspected in the same
way as described in the chapters \ref geometry_s2 and \ref geometry_s3. 

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

To be run from aliroot:
<pre>
AliMpCDB::LoadMpSegmentation2();
.x MUONCheckMisAligner.C
</pre>

It uses AliMUONGeometryAligner : 
- Creates a new AliMUONGeometryTransformer and AliMUONGeometryAligner
- Loads the geometry from the specified geometry file (default is geometry.root)
- Creates a second AliMUONGeometryTransformer by misaligning the existing 
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


\section geometry_s8 Geometry data files format
 
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


This chapter is defined in the READMEgeometry.txt file.

*/
