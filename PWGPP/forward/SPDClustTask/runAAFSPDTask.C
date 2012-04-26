//
void runAAFSPDTask(TString dataset="/alice/sim/LHC11d2_000119161",
		   TString outFName = "spdres.root",
		   Int_t   nEvents    = -1,//3000,
		   Float_t etaMin     =-1.5,        // min eta range to fill in histos
		   Float_t etaMax     = 1.5,        // max eta range to fill in histos
		   Float_t zMin       = -5,         // process events with Z vertex min
		   Float_t zMax       =  5,         //                     max positions
		   //
		   Float_t cutSigNStd  = 1.,       // cut on weighed distance used to extract signal
		   Float_t cutSigDPhiS = -1,        // cut on dPhi-phiBent used to extract signal (if negative -> dphi*sqrt(cutSigNStd)
		   Bool_t  useMC  = kTRUE,          // fill MC info (doRec=kTRUE)
		   //
		   // specific parameters for reconstruction
		   float  phiRot      = 3.14159e+00, // angle for bg. generation with rotation
		   Bool_t scaleDTheta = kTRUE,       // scale dTheta by 1/sin^2(theta) in trackleting
		   float  nStdDev     = 1.,         // number of st.dev. for tracklet cut to keep
		   float  dphi        = 0.06,        // dphi window (sigma of tracklet cut)
		   float  dtht        = 0.025,       // dtheta .... (if negative, abs will be used with additional cut on |dthetaX|, apart from w.distance
		   float  phishift    = 0.0045,      // bending shift
		   Bool_t remOvl      = kTRUE,       
		   float  ovlPhiCut   = 0.005, 
		   float  ovlZetaCut  = 0.05,
		   Int_t  nEventsSkip = 0,
		   //
		   TString alirootVer = "VO_ALICE@AliRoot::v5-03-15-AN",
		   TString rootVer    = "VO_ALICE@ROOT::v5-33-02a",
		   //
		   TString proofCluster="cheshkov@alice-caf.cern.ch"
		   ) 
{ 
  //  
  Bool_t runLocal = kFALSE;//kTRUE; // true only for local test mode
  if (runLocal) {
    //    dataset = "/default/shahoian/test_pp";//"/default/shahoian/test";
    dataset = "default/shahoian/test_sim_lhc11b1a";
    proofCluster = "";
    alirootVer = "AliRootProofLite";
    nEvents = 500;
  }
  //
  if ((!dataset.Contains("sim")) && useMC) {
    printf("Running with read data dataset, switching OFF useMC\n");
    useMC = kFALSE;
  }
  //
  printf("Requested: %s %s\n",alirootVer.Data(), rootVer.Data());
  printf("Output expected in %s\n",outFName.Data());
  //
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  //
  TString alirootMode="REC";//"REC";
  TString extraLibs = "ITSrec:CDB:Geom:"; // not needed in default aliroot mode
  //extraLibs+= "ANALYSIS:ANALYSISalice";
  extraLibs+= "ANALYSIS:OADB:ANALYSISalice:EventMixing";
  TList *list = new TList();
  // sets $ALIROOT_MODE on each worker to let proof to know to run in special mode
  list->Add(new TNamed("ALIROOT_MODE"      , alirootMode.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES", "ITS:include"));
  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN","1"));
  //
  //REM: same version of AliRoot on client!!!!! Otherwise error!! 
  TProof::Mgr(proofCluster.Data())->SetROOTVersion(rootVer.Data());
  //  TProof::Open(proofCluster.Data());//,"workers=10x");
  TProof::Open(proofCluster.Data(),"workers=1x");
  if (!gProof) {
    Error("runAAFMulti.C","Connection to AF failed.");
    return;
  }
  //  gProof->Exec("TObject *o = gEnv->GetTable()->FindObject(\"Proof.UseMergers\");"
  //	       "gEnv->GetTable()->Remove(o);", kTRUE);
  gProof->SetParameter("PROOF_UseMergers", 0);
  // Lets enable aliroot + extra libs on proof cluster
  if (runLocal) gProof->UploadPackage(alirootVer.Data());
  gProof->EnablePackage(alirootVer.Data(), list);
  //  gProof->EnablePackage(alirootVer.Data());
  //
  if (runLocal) {
    Int_t numWorkers = gProof->GetParallel();
    if (numWorkers<1) {printf("No workers\n"); return;}
    gProof->SetParameter("PROOF_PacketizerStrategy", (Int_t)0);
    int frac = (Int_t) 5 / numWorkers;
    if (frac<1) frac = 1;
    gProof->SetParameter("PROOF_PacketAsAFraction", frac);
  }
  //
  gROOT->LoadMacro(Form("%s/AnalysisSPDClustTask.C",gSystem->pwd()));
  TStopwatch sw;
  sw.Start();
  AnalysisSPDClustTask(dataset,outFName,nEvents,etaMin,etaMax,zMin,zMax,
		       cutSigNStd,cutSigDPhiS,useMC,
		       scaleDTheta,nStdDev,dphi,dtht,
		       phishift,remOvl,ovlPhiCut,ovlZetaCut,nEventsSkip);
  //
  sw.Stop();
  sw.Print();
}
