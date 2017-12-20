///
/// \file TestEMCALSimulation.C
/// \ingroup EMCAL_TestSimRec
/// \brief Simple macro to test EMCAL simulation
///
/// Simple macro to test EMCAL simulation.
/// It will take the Config.C macro that is sitting in the same place as the execution is performed.
/// In order to execute this you can do 
///    * Root5: aliroot -q -b -l $ALICE_ROOT/EMCAL/macros/TestSimuReco/TestEMCALSimulation.C
///    * Root6: aliroot -q -b -l $ALICE_ROOT/EMCAL/macros/TestSimuReco/LoadLibForConfig.C $ALICE_ROOT/EMCAL/macros/TestSimuReco/TestEMCALSimulation.C
///
/// Or directly in the root prompt 
///    root [1] .x LoadLibForConfig.C //Root6
///    root [2] .x TestEMCALSimulation.C
///
/// In order to find all the included classes in the Config.C one should add to the rootlogon.C file some paths
/// gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/  -I$ALICE_ROOT/include -I$ALICE_ROOT/ANALYSIS/macros -I$ALICE_ROOT/STEER -I$ALICE_ROOT/STEER/STEER -I$GEANT3DIR/include -I$GEANT3DIR/include/TGeant3");
/// or do it in the root prompt before execution.
///
/// \author : Jenn Klay, LLNL.
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS). 
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TStopwatch.h>
#include <TSystem.h>

#include "AliSimulation.h"
#include "AliLog.h"

#endif

#ifndef TRANSPORTMODEL
#define TRANSPORTMODEL

TString kTransportModel = "None"; // Set it in LoadLibForConfig.C

#endif //TRANSPORTMODEL

///
/// Main execution method
///
/// \param nev: number of events to be generated
/// \param raw: generate also raw data from digits?
///
void TestEMCALSimulation(Int_t nev =10, Bool_t raw = kFALSE)
{
  if(kTransportModel=="" || kTransportModel=="None")
  {
    AliLog::Message(AliLog::kInfo, 
                    "DO NOTHING *** Remember to load before LoadLibForConfig.C!! Set there the transport model!! ***", 
                    "TestEMCALSimulation.C", "TestEMCALSimulation.C", "TestEMCALSimulation()","TestEMCALSimulation.C", __LINE__);
    return;
  }
  
  AliSimulation simulator;
  simulator.SetConfigFile("Config.C");
  simulator.SetMakeSDigits("EMCAL");
  simulator.SetMakeDigits ("EMCAL");
  
  //simulator.SetRunGeneration(kFALSE); // Generate or not particles
  //simulator.SetRunSimulation(kFALSE); // Generate or not HITS (detector response) or not, start from SDigits
  
  if(raw)  simulator.SetWriteRawData("EMCAL","raw.root",kTRUE);
  
  //OCDB settings
  simulator.SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
                               Form("local://%s",gSystem->pwd()));
  
  // In case of anchoring MC, comment previous OCDB lines
  // select the appropriate year
  //simulator.SetDefaultStorage("alien://Folder=/alice/data/2011/OCDB");
  //simulator.UseVertexFromCDB();
  //simulator.UseMagFieldFromGRP();
  
  //Avoid the HLT to run
  simulator.SetRunHLT("");
  
  //Avoid QA
  //simulator.SetRunQA(":");
  
  TStopwatch timer;
  timer.Start();
  
  //  simulator.SetRunNumber(159582); // LHC11d run
  
  simulator.Run(nev);
  
  timer.Stop();
  timer.Print();
  
}
