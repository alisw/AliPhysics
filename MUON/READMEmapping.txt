// $Id$

/*! 

\page README_mapping MUON Mapping
 

See the detailed description of the mapping package in ALICE-INT-2003-025.
Since that time the mapping has been extended for slat and trigger
chamber segmentations an later on also to hold the description of the
top level connections of detection elements, including information
about DDLs, bus patches and also trigger configuration.

\section mapping_s0  Mapping in OCDB
  
The mapping is everywhere in the MUON code loaded from
OCDB. In case the mapping data files are changed, the mapping
in OCDB has to be regenerated in this way:

For changes in mapping/data:
<pre>
$> rm $ALICE_ROOT/OCDB/MUON/Calib/MappingData/Run0_999999999_v0_s0.root
$> cd $ALICE_ROOT/MUON
$> aliroot
root [0] AliMpCDB::WriteMpData(); 
</pre>

For changes in mapping/data_run:
<pre>
$> rm $ALICE_ROOT/OCDB/MUON/Calib/MappingRunData/Run0_999999999_v0_s0.root
$> cd $ALICE_ROOT/MUON
$> aliroot
root [0] AliMpCDB::WriteMpRunData(); 
</pre>

Note that mapping has to be loaded from OCDB almost each time
when using MUON classes; the loading of mapping depends on
the CBD manager state (the current run number, storage ...).
The standard way of loading mapping expects the CDB manager
in a state well defined beforehand; this way is used in
the MUON code:
<pre>
if ( ! AliMpCDB::LoadDDLStore() ) {
  AliFatal("Could not access mapping from OCDB !");
}
</pre>
In the same way, you can load AliMpManuStore, when manu serial 
numbers are needed:
<pre>
if ( ! AliMpCDB::LoadManuStore() ) {
  AliFatal("Could not access run-dependent mapping from OCDB !");
}
</pre>
In the interactive Root session, in case the CDB manager state is 
not defined, you can load mapping from the local OCDB files in this 
way:
<pre>
root [0] AliMpCDB::LoadDDLStore2();
root [2] AliMpCDB::LoadManuStore2();
</pre>

\section mapping_s1  Graphical User Interface
  
To use the GUI to plot DE segmentation run:

<pre> 
AliMpCDB::LoadDDLStore2();
AliMpCDB::LoadManuStore2();
new AliMpDEVisu();
</pre>

or

<pre>
AliMpCDB::LoadDDLStore2();
AliMpCDB::LoadManuStore2();
new AliMpDEVisu(w, h);
</pre>

if you want to change the size of the GUI window.
Typical value are:
<pre>
  w = 1200, h = 600 for PC
  w = 1000, h = 550 for laptop
</pre>

The GUI allows:
- drawing motif of a slat/quadrant
- search of a given manu (motif) number
- draw the channel number for a given manu number by clicking of the motif in canvas
- write down in log message informations about the given detection element
  * DE Id, DE name, 
  * number of buspatches, manus, manu serials
- option to save log message onto disc

\section mapping_s2 Test macros

A set of tests macros have been written during the development
of the mapping classes. To run these macros:

<pre>
   cd ../mapping/macro
   root
   root [0] .x testMacroName.C    
</pre>
                   

\section mapping_s3  Data files in mapping/data

The directory data in $ALICE_ROOT/MUON/mapping contains files
with data which are not supposed to be changed in a long period.

\subsection mapping_s3_sub1  zones.dat

Describes layout of zones, rows, row segments, subzones, motifs

<pre>
  SECTOR_DATA
    number of zones  
    number of rows  
    direction of constant pad size (X or Y)
    offset in X direction
    offset in Y direction
  
  ZONE     
    number of zone  
    half legth of pad size in x  
    half legth of pad size in y

  SUBZONE  
    motif id  
    motif type_id       

  ROW_SEGMENT  
    x offset (in number of pads) 
    y offset (in number of pads) 
    row number 
    nof motifs 
    first motif position Id
    step to the next motif position Id (+1 or -1)
</pre>
  

\subsection mapping_s3_sub2  zones_special.dat

Describes layout of special row segments (with irregular motifs)

<pre>
  SECTOR_SPECIAL_DATA

  MOTIF
    zone id
    motif id  
    motif type_id       

  ROW
    row number
  
  PAD_ROWS
    number of these pad rows in row   
  
  PAD_ROW_SEGMENT
    mumber of pads in the rows segment  
    motif id  
    motif position id
  
  motifX.dat
  ----------
  Describes characteristics of the motif type X

  In lines:
    Berg number
    Kapton number
    Pad number
    Gassi number
</pre>


\subsection mapping_s3_sub3  motifSpecialX.dat

Describes characteristics of the special motif with motif Id X;
the special motif caontains pads of different size

<pre>
  In lines:
    pad index i (in x)
    pad index j (in y)
    half legth of pad size in x  
    half legth of pad size in y
</pre>
  
\subsection mapping_s3_sub4  padPosX.dat

Maps pad numbers used in the motifX.dat files to
the local pad indices (i,j)

<pre>
  In lines:
    Pad number
    pad index i (in x)
    pad index j (in y)
</pre>
  

\subsection mapping_s3_sub5  *.pcb files

Lines starting with # are comments.

<pre>
  SIZES PadSizeX PadSizeY SizeX SizeY (cm)

  MOTIF motifType ix iy
  MOTIF motifType ix iy
  ...
</pre>

where ix, iy are the local coordinates (in pad unit) of the
lower-left corner of the motif (0,0 is the lower-left corner
of the PCB).

PCB *MUST* be described in a rotating way, starting lower-left and 
then counter-clockwise, otherwise the manu-to-motif association 
(fixed in the slat definition files) will be wrong.

Note that for "full" PCBs, the SizeX and SizeY are redundant as they could be 
computed from the motif alone (but that serves as a cross-check that the motif 
pattern given is ok). That's not the case for short or rounded PCB though.


\subsection mapping_s3_sub6  *.slat files

A slat is defined by the list of its PCB, described starting 
from the beam and going outward.

One PCB per line, preceded by the keyword PCB
Other lines not matching this syntax are ignored.
After the PCB is the list of manu ids for this PCB.

Example :

<pre>
  PCB X 1-3;24-20;42;44;53
  PCB X 1-14
  PCB Y 100-90
  PCB Z 1;2;3;4;5;6;7;12;120
</pre>

defines a slat with 4 PCBs : XXYZ

The manu to motif relationship is attached to the fact that we're counting 
counter-clockwise, starting on the lower-left of the PCB. (and the pcb files 
have to follow this convention to defined their motifs, otherwise all 
this won't work).

Note that the definition of the PCBs have to be in files with extension
.pcb (X.pcb, Y.pcb, Z.pcb)

  
\subsection mapping_s3_sub7  DetElemIdToBusPatch.dat

Lines starting with # are comments.

Contains the detection element identifier with the associated buspatch numbers 
and the corresponding DDL identifier.
The link between buspatches and DE's is needed on the rawdata level to identify 
the type of quadrant/slat to get the corresponding mapping.
The DDL id is needed for the rawdata generation only.

To generate this file, the macro MUONGenerateBusPatch.C could be used.


\subsection mapping_s3_sub8  BusPatchSpecial.dat

Lines starting with # are comments.

Contains the list of bus patches which manu readout is
not in the standard order. The format:

<pre>
KEYWORD  DDLs  BusPatches [ManuIDs}
         where KEYWORD = REVERT or EXPLICIT
</pre>

- For the bus patches following the REVERT keyword,
  the manus are just reordered in a reverted order.
- For the bus patches following the EXPLICIT keyword,
the manus filled with a standard procedure (using the DetElemIdToBusPatch.dat
file) are replaced with the list of manus in this file.


\subsection mapping_s3_sub9  BusPatchLength.dat

Lines starting with # are comments.

Contains the list of bus patches and their cable length in meters

<pre>
# DDL 0
  1  3
  2  3
...
</pre>


\subsection mapping_s3_sub10  crate.dat
  
Muon trigger electronics configuration file (decoded in class 
AliMUONTriggerCrateStore) directly copy/paste from the ALICE PRR 
ALICE-EN-2003-010. Gives local board number, name, 
crate name it belongs to, slot number, and internal switches 
(used in the algorithm).


\subsection mapping_s3_sub11  ManuSerialToBin.dat

Lines starting with # are comments.

Contains the manu serial number with their associated bin number, injection and calibration gain.

To compare the bin number with the serial in the CDB database you can run the macro:

<pre> 
AliMpCDB::LoadDDLStore2();
.L $ALICE_ROOT/MUON/mapping/macros/MUONCheckManu.C+
MUONCheckManu(10, kFALSE);
</pre>

The function has two parameters, the first is the number of the chamber (zero mean all chambers).
The macro can create a set of histogramms with the different gain distributions stored into a root file 
(second parameter).

Two files a generated: one with the list of manu per detection element with their associated bin and gain 
value, the other with the bad, strange or unidentified serial number.


\section mapping_s4  Data files in mapping/data_run

The directory data_run in $ALICE_ROOT/MUON/mapping contains files
with data which are expected to change during experiment.
At present time, there are only files with manu serial numbers:

\subsection mapping_s4_sub1  deName_manu_dat.dat

Contains the list of manuIds and their serial numbers.

<pre> 
# Id bp/nbp serial
1 bp 4615
2 bp 4616
...
</pre> 

\section mapping_s5  Units used
 
Lengths are in centimeters.
 
This chapter is defined in the READMEmapping.txt file.

*/
  
