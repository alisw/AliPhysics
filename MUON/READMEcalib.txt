// $Id$

/*! 

\page README_calib Calibration
 
The Offline Condition DataBase is described extensively on ALICE Offline pages.

Here you'll find only information relevant to the AliMUONCDB class
 (formerly MUONCDB.C macro), which defines a set of functions to read/write 
 MUON information to this CDB. Those functions are not meant to be used 
 as black boxes.

Please have a closer look before using (especially the ones writing to the CDB...)
 

\section calib_s1 Calibration data objects

We've designed generic data containers to store calibration information, 
tailored to the way we usually access MUON tracker data, that is, 
indexed by the pair (detElemId,manuId). 
This container is called AliMUONV2DStore. You can attach a TObject to every and
each pair (detElemId,manuId).

For the moment, that TObject is generally of AliMUONVCalibParam type, 
 which handles a given number of channels (64 typically) as a group. 
There's also an AliMUONV1DStore for data types which only requires indexing 
by 1 value (like trigger masks for instance).

As the class names suggest (V...), those classes are only interfaces. 
Concrete ones are AliMUON2DMap (used instead of a vector as detElemId are 
not contiguous) for the V2DStore, AliMUON1DArray (for things where indices are
contiguous) and AliMUON1DMap for the V1DStore, and CalibParamNI (VCalibParam 
storing n integer per channel), and CalibParamNF 
(VCalibParam storing n floats per channel).

One exception are the HV values from DCS, which are stored "as they come" from 
the shuttle-dcs interface, as a TMap, where the key is the aliasname (TString), 
and the value a TObjArray of AliDCSValue.

For trigger, the same virtual container idea applies, 
except we're using 1D container (AliMUONV1DStore, for masks) or specific ones (for LUT
and efficiency)

\section calib_s2 CDB location

One very important notion is that of the DefaultStorage (which you might set with 
 AliCDBManager::Instance()->SetDefaultStorage(path)), which tells the CDB library where
 the CDB is sitting (either locally on disk, or on the grid).

For local tests, path will be most likely = <code> "local://$ALICE_ROOT/OCDB"</code>
(i.e. there is, in CVS a slim version of the calibration objects needed
for running the MUON code), or <code> "local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestCDB </code>
for Shuttle testing.

When using alien, the path definition can be eg.
<code> "alien://folder=/alice/data/2007/LHC07w/OCDB" </code>. 
 
\section calib_s3 Writing to CDB
 
AliMUONCDB class is used to populate the CDB with fake calibration objects for testing purposes.
Real calibration data will normally be handled by the Shuttle (see READMEshuttle).

The various WriteXXX() methods may be used to populate MUON tracker and trigger 
information.
A full set of calibration data types can be created from scratch using, from
the Root prompt (from within the $ALICE_ROOT/MUON directory to get the correct
list of libraries loaded by the loadlibs.C macro)

<pre>
root[0] AliMpCDB::LoadDDLStore2(); 
root[1] AliMpCDB::LoadManuStore2(); 
root[2] AliMUONCDB cdb;
root[3] Int_t startRun = 0;
root[4] Bool_t defaultValues = kTRUE;
root[5] cdb.WriteTrigger(startRun);
root[6] cdb.WriteTracker(defaultValues,startRun);
</pre>


\section calib_s4 Reading the CDB
 
The actual reading is encapsulated into AliMUONCalibrationData class. 
e.g. to read pedestals for run 4567, one would do :

<pre>
AliCDBManager::Instance()->SetDefaultStorage(cdbPath);
AliMUONCalibrationData cd(4567);
AliMUONV2DStore* ped = cd.Pedestals();
</pre>

If you want to plot calibration data (not terribly usefull as it's a really global view),
 use the Plot() function of AliMUONCDB, e.g.  
 
<pre>
AliMUONCDB cdb(cdbpath);
cdb.Plot(*ped,"pedestal")
</pre>

which will create 2 histograms : pedestal_0 (mean) and pedestal_1 (sigma).

You might also be interested in the AliMUONStoreHelper::Diff() method 
which generates an AliMUONV2DStore containing the difference 
(either absolute or relative) of two AliMUONV2DStore.

\section calib_s5 A note on status and status map

A special kind of object is the status map. It would deserve a full documentation
 (and that will need to be done some day), but for the moment, please have a look
 at the MUONStatusMap.C macro which can let you play with it. For more information,
 have a look at the AliMUONPadStatusMaker and AliMUONPadStatusMapMaker.

This chapter is defined in the READMEcalib.txt file.
*/

