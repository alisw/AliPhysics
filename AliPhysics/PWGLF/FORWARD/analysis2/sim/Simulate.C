/**
 * @file   Simulate.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Oct 15 13:28:09 2014
 * 
 * @brief  Steering script for the simulation 
 */
Bool_t RunLego(AliSimulation& steer)
{
  TString runType = gSystem->Getenv("CONFIG_RUN_TYPE");
  runType.ToUpper();
  if (!runType.BeginsWith("LEGO")) return false;

  TString v(runType);
  v.ToUpper();
  v.ReplaceAll("LEGO", "");
  Info("Setup", "Making Lego event generator (variant: %s)", v.Data());

  AliLegoGenerator* gener = 0;
  if (v.BeginsWith("X") || v.BeginsWith("Y") || v.BeginsWith("Z")) {
      const char* c[] = { v(0), '\0' };
      gener = new AliLegoGeneratorXYZ(c);
  }
  else if (v.BeginsWith("PHIZ")) {
    gener = new AliLegoGeneratorPhiZ();
  }
  else if (v.BeginsWith("ETA")) {
    gener = new AliLegoGeneratorEta();
  }
  else {
    gener = new AliLegoGenerator();
  }
  
  // XYZ varies origin of the particles in two dimensions:
  //  X:  o=(0,t1,t2), p=(1,0,0)
  //  Y:  o=(t1,0,t2), p=(0,1,0)
  //  Z:  o=(t1,t2,0), p=(0,0,1)
  // PhiZ varies the momentum in two dimensions
  //  o=(0,0,t1) p=(cos(t2),sin(t2),0)
  // Eta varies momentum in two dimensions
  //  phi=t1
  //  theta=2*atan(exp(-t2))
  //  o=(0,0,0) p=(cos(phi)*sin(theta),sin(phi)*cos(theta),cos(theta))
  // Base varies in two dimensions
  //  phi=t1
  //  theta=t2
  //  o=(0,0,0) p=(cos(phi)*sin(theta),sin(phi)*cos(theta),cos(theta))
  const char* cfg = "Config.C";
  Bool_t ret = false;
  Double_t rMin = 0;
  Double_t rMax = 32; // 430;
  Double_t zMax = 400; // 10000;
  if (v.BeginsWith("X") || v.BeginsWith("Y"))
    ret = steer.RunLego(cfg,
			 10,-2,2,      // Y, X
			 10,-10,10,    // Z
			 rMin, rMax, zMax, gener);
  else if (v.BeginsWith("Z"))
    ret = steer.RunLego(cfg,
			 10,-2,2,      // X
			 10,-2,2,      // Y
 			 rMin, rMax, zMax, gener);
  else if (v.BeginsWith("PHIZ"))
    ret = steer.RunLego(cfg,
			 10,-10,10, // Z
			 360,0,360, // phi
			 rMin, rMax, zMax, gener);
  else if (v.BeginsWith("ETA")) {
    Double_t aEta = 6;
    Double_t dEta = (6.--4.)/200;
    ret = steer.RunLego(cfg,
			360,0,360, // phi
			2*aEta/dEta,-aEta,+aEta, // Eta
			rMin, rMax, zMax, gener);
  }
  else 
    ret = steer.RunLego(cfg);
  if (!ret)
    Warning("RunLego", "Failed to do lego run");
  return true;
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
  if (detCfg->GeometrySource()) 
    steer.SetGeometryFile(detCfg->GeometrySource());

  // -----------------------------------------------------------------
  //
  // OCDB and specific storages 
  // 
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetRun(grp->run);
  cdb->SetDefaultStorageFromRun(grp->run);
  cdb->SetRun(-1);
  ocdbCfg->Init(true);
  steer.SetRunNumber(grp->run);
  
  // -----------------------------------------------------------------
  // 
  // The rest - disable QA and HLT (memory heavy) for PbPb
  //
  if (isAA) steer.SetRunQA(":");
  if (is10h) steer.SetRunHLT("");
  
  TStopwatch timer;
  timer.Start();
  // Check first if we're doing a Lego run, and if not,
  // do a normal run
  if (!RunLego(steer)) steer.Run(nev);
  timer.Stop();
  timer.Print();
}
// 
// EOF
//  
