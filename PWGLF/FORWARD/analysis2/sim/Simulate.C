/** 
 * Read an integer valued environment variable 
 * 
 * @param name environment variable name 
 * 
 * @return Value or 0
 */
Int_t getIntEnv(const char* name)
{
  TString env = gSystem->Getenv(name);
  if (env.IsNull()) return 0;
  return env.Atoi();
}
/** 
 * Set specific storage 
 * 
 * @param s   Simlation object
 * @param key Key to set specific storage for
 * @param sub Sub-component of storage
 */
void SetSpecStore(AliSimulation& s, 
		  const char* key, 
		  const char* sub)
{
  s.SetSpecificStorage(key, Form("alien://Folder=/alice/simulation/%s",sub));
}

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
  // Get GRP parameters.  Defines global "grp" as a pointer to GRPData
  //
  gROOT->Macro(Form("GRP.C(%d)", run));

  // --- Get GRP to deduce collision system --------------------------
  Bool_t         isAA  = grp->IsAA();
  Bool_t         isPP  = grp->IsPP();
  Bool_t         is10h = grp->period.EqualTo("LHC10h");

  // -----------------------------------------------------------------
  // 
  // Basic setup 
  //
  AliSimulation steer; 
  steer.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  steer.SetMakeDigitsFromHits("ITS TPC");

  // -----------------------------------------------------------------
  // 
  // Vertex, Mag.field, and trigger from OCDB
  //
  steer.SetTriggerConfig(!isAA ? "p-p" : "Pb-Pb");//Replace with "ocdb"
  steer.UseMagFieldFromGRP();
  steer.UseVertexFromCDB();

  // -----------------------------------------------------------------
  //
  // Raw OCDB
  // 
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorageFromRun(grp->run);
  // cdb->SetRun(grp.run);
  steer.SetDefaultStorage(cdb->GetDefaultStorage()->GetURI());

  // --- ITS  (1 Total) ----------------------------------------------
  SetSpecStore(steer,"ITS/Align/Data",	"2008/v4-15-Release/Ideal");
  
  // --- MUON (1 object) ---------------------------------------------
  SetSpecStore(steer,"MUON/Align/Data",	"2008/v4-15-Release/Ideal"); 

  // ---- TPC (6 total) ----------------------------------------------
  SetSpecStore(steer,"TPC/Calib/TimeGain",	"2008/v4-15-Release/Ideal/");
  SetSpecStore(steer,"TPC/Calib/ClusterParam",	"2008/v4-15-Release/Ideal/");
  SetSpecStore(steer,"TPC/Calib/AltroConfig",	"2008/v4-15-Release/Ideal/");
  SetSpecStore(steer,"TPC/Calib/Correction",	"2008/v4-15-Release/Ideal/");
  SetSpecStore(steer,"TPC/Align/Data",		"2008/v4-15-Release/Ideal/");
  SetSpecStore(steer,"TPC/Calib/RecoParam",	"2008/v4-15-Release/Residual");
  if (is10h)
    SetSpecStore(steer,"TPC/Calib/TimeDrift",	"2008/v4-15-Release/Residual/");
  else 
    SetSpecStore(steer,"TPC/Calib/TimeDrift",	"2008/v4-15-Release/Ideal/");
    
  // --- ZDC for 2010 the following is needed ------------------------
  // (https://savannah.cern.ch/task/?func=detailitem&item_id=33180#comment46)
  if (is10h)
    SetSpecStore(steer,"ZDC/Align/Data",	"2008/v4-15-Release/Ideal/"); 


  // -----------------------------------------------------------------
  // 
  // The rest - disable QA for PbPb
  //
  if (isAA) steer.SetRunQA(":");
  // gInterpreter->UnloadFile("GetGRP.C");

  TStopwatch timer;
  timer.Start();
  steer.Run(nev);
  timer.Stop();
  timer.Print();
}
// 
// EOF
//  
