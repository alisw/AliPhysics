/**
 * @file   Simulate.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Oct 15 13:28:09 2014
 * 
 * @brief  Steering script for the simulation 
 */
/** 
 * Run the simulation 
 * 
 * @param nev Number of events per job
 * @param run Run number to simulate 
 */
void Simulate(Int_t nev=1, UInt_t run=0) 
{
  // -----------------------------------------------------------------
  // 
  // - Get GRP parameters.  Defines global "grp" as a pointer to GRPData
  // - Load base class definitions in BaseConfig.C
  // - Get which detectors are turned on in "detCfg". 
  // - Create the OCDB configuration object "ocdbCfg"
  // 
  gROOT->Macro(Form("GRP.C(%d)", run));
  gROOT->Macro("BaseConfig.C");
  gROOT->Macro("DetConfig.C"); 
  gROOT->Macro("OCDBConfig.C"); 

  // --- Get GRP to deduce collision system --------------------------
  Bool_t         isAA  = grp->IsAA();
  Bool_t         isPP  = grp->IsPP();
  Bool_t         is10h = grp->period.EqualTo("LHC10h");

  // -----------------------------------------------------------------
  // 
  // Basic setup 
  //
  AliSimulation steer; 
  TString sDigits, fromHits;
  detCfg->GetSDigitString(sDigits);
  detCfg->GetHits2DigitsString(fromHits);
  steer.SetMakeSDigits(sDigits);
  steer.SetMakeDigitsFromHits(fromHits);

  // -----------------------------------------------------------------
  // 
  // Vertex, Mag.field, and trigger from OCDB
  //
  steer.SetTriggerConfig(!isAA ? "p-p" : "Pb-Pb");//Replace with "ocdb"
  steer.UseMagFieldFromGRP();
  steer.UseVertexFromCDB();

  // -----------------------------------------------------------------
  //
  // OCDB and specific storages 
  // 
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorageFromRun(grp->run);
  ocdbCfg->Init(true);

  // -----------------------------------------------------------------
  // 
  // The rest - disable QA and HLT (memory heavy) for PbPb
  //
  if (isAA) steer.SetRunQA(":");
  if (is10h) steer.SetRunHLT("");
  
  TStopwatch timer;
  timer.Start();
  steer.Run(nev);
  timer.Stop();
  timer.Print();
}
// 
// EOF
//  
