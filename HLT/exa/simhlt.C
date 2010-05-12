// $Id$
/**
 * @file simhlt.C
 * @brief Run HLT reconstruction embedded into AliSimulation
 *
 * <pre>
 * Usage: aliroot -b -q -l \
 *     aliroot -b -q -l simhlt.C'("hltoption", "rawdataoptions", nofEvents, runNo, "cdbUri")'
 *
 * Examples:
 *     simhlt.C'("loglevel=0x7c chains=GLOBAL-esd-converter")'
 *       load all default libraries and run the chain 'GLOBAL-esd-converter'
 *     simhlt.C'("loglevel=0x7c !libAliHLTITS.so")'
 *       load all default libraries but the libAliHLTITS.so and run the
 *       default chain, skips the ITS tracking
 *     simhlt.C'("loglevel=0x7c rawfile=")'
 *       run explicitly on digit data
 *
 * Defaults
 *     rawdataoptions=""         -> skip generation of raw data
 *     nofEvents=-1              -> take event count from simulated sample
 *     runNo=-1                  -> take run number from simulated sample
 *     cdbUri="local://$ALICE_ROOT/OCDB"
 *
 * </pre>
 *
 * This macro can be used to run the 'HLT simulation' standalone with an
 * already simulated data sample. 'HLT simulation' means HLT reconstruction
 * on simulated data, and is carried out at the end of the AliSimulation
 * sequence.
 *
 * The definition of the HLT reconstruction chain depends on the available
 * module plugins (component libraries). Each library defines a number of
 * component configurations identified by a unique string. Here are a couple
 * of definitions, but the list is by far not complete:
 * - libAliHLTTPC.so:     TPC-clusters, TPC-globalmerger
 * - libAliHLTITS.so:     ITS-tracker, ITS-SPD-vertexer
 * - libAliHLTGlobal.so:  GLOBAL-esd-converter, GLOBAL-vertexhisto
 * - libAliHLTTrigger.so: GLOBAL-Trigger
 *
 * This list will certainly develop over time, it is planned to add some
 * helper functionality to AliHLTSystem for querying of such information.
 *
 * The definition of the chain depends on the availability of simulated raw
 * data. If raw data is available, the HLT chain runs on this, otherwise on
 * the digit data. \b NOTE: propagation of MC information is only possible
 * for the latter case. Raw data generation must be either switched of, i.e.
 * rawdataoptions is empty, or explicitly explicitly ignored by the HLT
 * option 'rawfile=' (note the empty argument to the option).
 *
 * The output of the simulation are the HLT.Digits.root file and (if enabled)
 * the HLT raw data files. You must run AliReconstruction in order to extract
 * the simulated HLT data. A typical test sequence can be like:
 *
 * <pre>
 * cd $ALICE_ROOT/test/ppbench
 * ./runtest.sh
 * aliroot -b -q -l simhlt.C'("loglevel=0x7c chains=GLOBAL-esd-converter")' \
 *         | tee simhlt.C
 * aliroot -b -q -l rec.C | tee rec.log
 * </pre>
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void simhlt(const char *hltOptions,
	    const char *rawDataSelection="",
	    int nofEvents=-1,
	    int runNo=-1,
	    const char *cdbURI="local://$ALICE_ROOT/OCDB"
	    )
{
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // connect to the GRID if we use a file or OCDB from the GRID
  TString struri=cdbURI;
  if (struri.Contains("://") && !struri.Contains("local://")) {
    TGrid::Connect("alien");
  }

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbURI);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the HLT simulation
  // All but HLT simulation is switched off
  //
  AliSimulation sim;
 
  // switch of simulation and data generation
  // comment all that stuff to also simulate the events and data
  sim.SetRunGeneration(kFALSE);
  sim.SetMakeDigits("");
  sim.SetMakeSDigits("");
  sim.SetMakeDigitsFromHits("");
  //sim.SetMakeTrigger("");
  sim.SetRunQA("HLT:ALL");

  // the normal simulation sets the specific storage for the GRP entry
  if (gSystem->AccessPathName("GRP/GRP/Data")) {
    cerr << "*********************************************************" << endl;
    cerr << "error: no GRP entry found in the currect directory, simulation might be incomplete. Skip setting specific storage for GRP entry" << endl;
    cerr << "*********************************************************" << endl << endl;
  } else {
    sim.SetSpecificStorage("GRP/GRP/Data",
			   Form("local://%s",gSystem->pwd()));
  }

  // set the options for writing raw data, possible values e.g. ALL; HLT
  if (rawDataSelection && rawDataSelection[0]!=0) {
    sim.SetWriteRawData(rawDataSelection, "raw.root", kTRUE);
  }

  // set the options for the HLT simulation
  sim.SetRunHLT(hltOptions);

  if (runNo>=0) sim.SetRunNumber(runNo);

  if (nofEvents>0) sim.Run(nofEvents);
  else sim.Run();
}

void simhlt()
{
  cout << "recraw-local: Run HLT reconstruction embedded into AliSimulation" << endl;
  cout << " Usage: aliroot -b -q -l \\" << endl;
  cout << "     aliroot -q simhlt.C'(\"hltoption\", \"rawdataoptions\", nofEvents, runNo, \"cdbUri\")'" << endl;
  cout << "" << endl;
  cout << " Examples:" << endl;
  cout << "     aliroot -q simhlt.C'(\"loglevel=0x7c chains=GLOBAL-esd-converter\")'" << endl;
  cout << "     aliroot -q simhlt.C'(\"loglevel=0x7c !libAliHLTITS.so\")'		 " << endl;
  cout << "     aliroot -q simhlt.C'(\"loglevel=0x7c rawfile=\")'                   " << endl;
  cout << "" << endl;
  cout << " Defaults" << endl;
  cout << "     rawdataoptions=\"\"         -> skip generation of raw data	   " << endl;
  cout << "     nofEvents=-1              -> take event count from simulated sample" << endl;
  cout << "     runNo=-1                  -> take run number from simulated sample " << endl;
  cout << "     cdbUri=\"local://$ALICE_ROOT/OCDB\"                                " << endl;
  cout << "" << endl;
  cout << " The definition of the chain depends on the availability of simulated raw " << endl;
  cout << " data. If raw data is available, the HLT chain runs on this, otherwise on " << endl;
  cout << " the digit data. \b NOTE: propagation of MC information is only possible  " << endl;
  cout << " for the latter case. Raw data generation must be either switched of, i.e." << endl;
  cout << " rawdataoptions is empty, or explicitly explicitly ignored by the HLT	  " << endl;
  cout << " option 'rawfile=' (note the empty argument to the option).               " << endl;
}
