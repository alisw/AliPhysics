// $Id$

/*! 

\page README_calib MUON Calibration
 
 
The Offline Condition DataBase is described extensively on ALICE Offline pages.

Here you'll find mainly information relevant to the MUON objects in the OCDB,
and to the AliMUONCDB class (formerly MUONCDB.C macro), which defines a set of functions to read/write 
 MUON information to this CDB. Those functions are not meant to be used 
 as black boxes.

Please have a closer look before using (especially the ones writing to the CDB...)

\section calib_s0 MUON objects in OCDB
 
<table>
<tr><td>OCDB path</td><td>Subsystem</td><td>Run Type</td><td>FileID</td><td>DAQsource</td><td>Update frequency</td><td>Appox. Size</td></tr>
<tr><td>MUON/Calib/RecoParam</td><td>both</td><td>N/A</td><td>N/A</td><td>N/A</td><td>Few per year ?</td><td>5 KB</td></tr>
<tr><td>MUON/Calib/MappingData</td><td>both</td><td>N/A</td><td>N/A</td><td>N/A</td><td>Once for all once debugged</td><td>1 MB</td></tr>
<tr><td>MUON/Calib/MappingRunData</td><td>both</td><td>N/A</td><td>N/A</td><td>N/A</td><td>Between zero and few per year</td><td>100 KB</td></tr>
<tr><td>MUON/Calib/HV</td><td>MCH</td><td>PHYSICS</td><td>N/A</td><td>N/A</td><td>Once per physics run</td><td>Depends on the number of trips, 10-20 KB normally</td></tr>
<tr><td>MUON/Calib/Neighbours</td><td>MCH</td><td>N/A</td><td>N/A</td><td>N/A</td><td>As MappingData</td><td>10 MB</td></tr>
<tr><td>MUON/Calib/OccupancyMap</td><td>MCH</td><td>PHYSICS</td><td>OCCUPANCY</td><td>MON</td><td>Once per physics run</td><td>Depends on the number of bad manus, normally 100-200 KB</td></tr>
<tr><td>MUON/Calib/OccupancyMap</td><td>MCH</td><td>PHYSICS</td><td>OCCUPANCY</td><td>MON</td><td>Once per physics run</td><td>Depends on the number of run duration and the time resolution used</td></tr>
<tr><td>MUON/Calib/Pedestals</td><td>MCH</td><td>PEDESTAL</td><td>PEDESTALS</td><td>LDC</td><td>Once per pedestal run</td><td>7 MB</td></tr>
<tr><td>MUON/Calib/RejectList</td><td>MCH</td><td>N/A</td><td>N/A</td><td>N/A</td><td>Few per year</td><td>Depends on the number of bad elements</td></tr>

<tr><td>MUON/Calib/GlobalTriggerCrateConfig</td><td>MTR</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td></tr>
<tr><td>MUON/Calib/LocalTriggerBoardMasks</td><td>MTR</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td></tr>
<tr><td>MUON/Calib/RegionalTriggerConfig</td><td>MTR</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td></tr>
<tr><td>MUON/Calib/TriggerDCS</td><td>MTR</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td></tr>
<tr><td>MUON/Calib/TriggerEfficiency</td><td>MTR</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td></tr>
<tr><td>MUON/Calib/TriggerLut</td><td>MTR</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td></tr>
</table>

In addition, the following ones were used in the past, but have been discontinued (but are still present in the OCDB, as object removal is not allowed).

<table>
<tr><td>OCDB path</td></tr>
<tr><td>MUON/Calib/DDLStore</td></tr>
<tr><td>MUON/CalibGlobalTriggerBoardMasks/</td></tr>
<tr><td>MUON/Calib/Mapping</td></tr>
<tr><td>MUON/Calib/RegionalTriggerBoardMasks</td></tr>
<tr><td>MUON/Calib/RegionalTriggerConfig?</td></tr>
</table>

\section calib_s1 Calibration data object classes

We've designed generic data containers to store calibration information, 
tailored to the way we usually access MUON tracker data, that is, 
indexed by the pair (detElemId,manuId). 
This container is called AliMUONVStore. You can attach a TObject to every and
each pair (detElemId,manuId).

For the moment, that TObject is generally of AliMUONVCalibParam type, 
 which handles a given number of channels (64 typically) as a group. 
There's also 1D versions of AliMUONVStore for data types which only requires indexing 
by 1 value (like trigger masks for instance).

As the class names suggest (V...), those classes are only interfaces. 
Concrete ones are AliMUON2DMap (used instead of a vector as detElemId are 
not contiguous), AliMUON1DArray (for things where indices are
contiguous) and AliMUON1DMap, and CalibParamNI (VCalibParam 
storing n integer per channel), CalibParamNF 
(VCalibParam storing n floats per channel), and CalibParamND (VCalibParam storing n doubles per channel).

One exception are the HV values from DCS, which are stored "as they come" from 
the shuttle-dcs interface, as a TMap, where the key is the aliasname (TString), 
and the value a TObjArray of AliDCSValue.

\section calib_s2 CDB location

One very important notion is that of the DefaultStorage (which you might set with 
 AliCDBManager::Instance()->SetDefaultStorage(path)), which tells the CDB library where
 the CDB is sitting (either locally on disk, or on the grid).

For local tests, path will be most likely = <code> "local://$ALICE_ROOT/OCDB"</code>
(i.e. there is, in svn a slim version of the calibration objects needed
for running the MUON code), or <code> "local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestCDB </code>
for Shuttle testing.

When using alien, the path definition can be eg.
<code> "alien://folder=/alice/data/2009/OCDB" </code>. 
 
\section calib_s3 Writing to CDB
 
AliMUONCDB class is used to populate the CDB with fake calibration objects for testing purposes.
Real calibration data will normally be handled by the Shuttle (see READMEshuttle).

The various WriteXXX() static methods may be used to populate MUON tracker and trigger 
information.
A full set of calibration data types can be created from scratch using, from
the Root prompt (from within the $ALICE_ROOT/MUON directory to get the correct
list of libraries loaded by the loadlibs.C macro)

<pre>
root[0] AliMpCDB::LoadAll2(); 
root[1] Int_t startRun = 0;
root[2] Bool_t defaultValues = kTRUE;
root[3] AliMUONCDB::WriteTrigger(startRun);
root[4] AliMUONCDB::WriteTracker(defaultValues,startRun);
</pre>


\section calib_s4 Reading the CDB
 
The actual reading is encapsulated into AliMUONCalibrationData class. 
e.g. to read pedestals for run 67138, one would do :

<pre>
AliCDBManager::Instance()->SetDefaultStorage(cdbPath);
AliMUONVStore* ped = AliMUONCalibrationData::CreatePedestals(67138);
</pre>

If you want to plot calibration data (not terribly usefull as it's a really global view),
 use the Plot() function of AliMUONCDB (this require to load the segmentation), e.g.  
 
<pre>
AliCDBManager::Instance()->SetRun(0);
AliMUONCDB::LoadMapping(kTRUE);
AliMUONCDB::Plot(*ped,"pedestal");
</pre>

which will create 2 histograms : pedestal_0 (mean) and pedestal_1 (sigma).

A more usefull way to look at calibration data might be the \link README_mchview mchview \endlink program.

\section calib_s5 A note on status and status map

A special kind of object is the status map. It would deserve a full documentation
 (and that will need to be done some day), but for the moment, please have a look
 at the MUONStatusMap.C macro which can let you play with it. For more information,
 have a look at the AliMUONPadStatusMaker and AliMUONPadStatusMapMaker.

 $Id$

*/

