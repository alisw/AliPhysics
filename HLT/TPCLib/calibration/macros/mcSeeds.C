// $Id$
/*
 * @file mcSeeds.C
 * @brief Calculating AliTPCseed objects for estimating dEdx.
 *
 * Example macro to run the HLT tracker embedded into
 * AliRoot simulation. The reconstruction is done from the TPC digits.
 * The seed making process carries the MC label information for comparison
 * of the HLT seeds with the offline ones.
 *
 * Usage: aliroot -b -q mcSeeds.C | tee seeds.log
 *
 * The macro assumes the data to be already simulated. If it should run
 * within the initial simulation, comment the corresponding functions
 * below (SetRunGeneration etc.)
 *
 * @author Kalliopi.Kanaki@ift.uib.no
 * @ingroup alihlt_tpc
 */

void mcSeeds(){
   
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //

  AliHLTSystem *gHLT = AliHLTPluginBase::GetInstance();

  // load TPCParam from OCDB
  AliCDBManager* pMan=AliCDBManager::Instance();
  if (pMan) {
    if (!pMan->IsDefaultStorageSet()) {
      pMan->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
      pMan->SetRun(0);
    }
  }
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
 
  TString seedMakerInput;
  if(seedMakerInput.Length()>0) seedMakerInput+=" ";
  seedMakerInput+="TPC-clusters TPC-globalmerger TPC-mcTrackMarker";

  AliHLTConfiguration seedconf("seeds", "TPCCalibSeedMaker", seedMakerInput.Data(), "");
  //AliHLTConfiguration dedxconf("dedx", "TPCdEdx", seedMakerInput.Data(), "");

  AliHLTConfiguration rfwconf("rfw", "ROOTFileWriter", "seeds", "-datafile hlt_seeds -idfmt=_%d");
  //AliHLTConfiguration rfwconf("rfw", "ROOTFileWriter", "seeds", "-datafile dedx.root -concatenate-events -overwrite");
  
  

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the HLT simulation
  // All but HLT simulation is switched off
  //
  AliSimulation sim;

  // switch off simulation and data generation
  // comment all that stuff to also simulate the events and data
  sim.SetRunGeneration(kFALSE);
  sim.SetMakeDigits("");
  sim.SetMakeSDigits("");
  sim.SetMakeDigitsFromHits("");
  //sim.SetMakeTrigger("");
  sim.SetRunQA(":");

  // the normal simulation sets the specific storage for the GRP entry
  if (gSystem->AccessPathName("GRP/GRP/Data")) {
    cerr << "*********************************************************" << endl;
    cerr << "error: no GRP entry found in the currect directory, simulation might be incomplete. Skip setting specific storage for GRP entry" << endl;
    cerr << "*********************************************************" << endl << endl;
  } else {
    sim.SetSpecificStorage("GRP/GRP/Data",  Form("local://%s",gSystem->pwd()));
  }

  // set the options for the HLT simulation
  sim.SetRunHLT("libAliHLTUtil.so libAliHLTTPC.so libAliHLTGlobal.so libAliHLTTPCCalibration.so loglevel=0x7c chains=rfw");
  sim.Run();

}



