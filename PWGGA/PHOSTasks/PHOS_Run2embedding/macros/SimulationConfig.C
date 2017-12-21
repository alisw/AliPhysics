/*
 * AliDPG - ALICE Experiment Data Preparation Group
 * Simulation configuration script
 *
 */

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

enum ESimulation_t {
  kSimulationDefault,
  kSimulationMuon,
  kSimulationCustom,
  kNSimulations
};

const Char_t *SimulationName[kNSimulations] = {
  "Default",
  "Muon",
  "Custom"
};


/*****************************************************************/

SimulationConfig(AliSimulation &sim, ESimulation_t tag, Int_t run)
{
  
  
  switch(tag) {
    
    // Default
  case kSimulationDefault:
    SimulationDefault(sim, run);
    return;
    
    // Muon
  case kSimulationMuon:
    SimulationDefault(sim, run);
    sim.SetMakeSDigits("MUON VZERO");
    sim.SetMakeDigitsFromHits("ITS");
    return;
    
    // Custom
  case kSimulationCustom:
    if ((gROOT->LoadMacro("SimulationCustom.C")) != 0) {
      printf("ERROR: cannot find SimulationCustom.C\n");
      abort();
      return;
    }
    SimulationCustom(sim, run);
    return;

  }
  
}

SimulationDefault(AliSimulation &sim, Int_t run)
{

  gROOT->LoadMacro("Sim/Utils.C");
  Int_t year = RunToYear(run);
  
  //
  // set OCDB snapshot mode
  //    AliCDBManager *man = AliCDBManager::Instance();
  //    man->SetDefaultStorage("alien://Folder=/alice/data/2015/OCDB");
  //    man->SetRun(run);
  //    man->SetSnapshotMode("OCDBsim.root");
  sim.SetCDBSnapshotMode("Sim/OCDBsim.root");
  //
//   if (year < 2015) sim.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON ZDC PMD T0 VZERO FMD");
//   else             sim.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON ZDC PMD T0 VZERO FMD AD");
  sim.SetMakeSDigits("PHOS");
  sim.SetMakeDigits("PHOS");
  
//  sim.SetMakeDigitsFromHits("ITS TPC");


  //  sim.SetRunHLT(""); // can't detect from GRP if HLT was running, off for safety now 
  //

	//HLT settings //this is copied from v5-08-XX-12, $ALIDPG_ROOT/SimulationCongig.C
	TString hltConfig = "auto";
	if (gSystem->Getenv("CONFIG_HLT"))
					hltConfig = gSystem->Getenv("CONFIG_HLT");
	sim.SetRunHLT(hltConfig.Data());

	//

  SimulationConfigPHOS(sim, run);
  //
  sim.UseVertexFromCDB();
  sim.UseMagFieldFromGRP();
  //
  sim.SetRunQA(":");
  //
}

/*** PHOS ****************************************************/

SimulationConfigPHOS(AliSimulation &sim, Int_t run)
{
  AliPHOSSimParam *simParam = AliPHOSSimParam::GetInstance();
	simParam->SetCellNonLinearity(kFALSE);//switch OFF cell non-linearity
  //simParam->SetCellNonLineairyA(-0.03);//from LHC15n
  //simParam->SetCellNonLineairyB(0.6);//from LHC15n
  //simParam->SetCellNonLineairyC(1.006);//from LHC15n
}

