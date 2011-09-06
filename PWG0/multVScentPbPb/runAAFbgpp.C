//
void runAAFbgpp(TString dataset="/alice/sim/LHC11d3_000146806", //"/alice/sim/LHC10f8c_130844",
		TString outFName = "LHC11d3_000146806_v1.root",
		Int_t   nEvents    = -1,//3000,
		TString noMergeDir = "root://alicers01.cern.ch//tmp/myoutput/", // "" // if non-zero string, no merging is done
		Bool_t scaleDTheta = kTRUE,       // scale dTheta by 1/sin^2(theta) in trackleting
		float  nStdDev     = 25.,         // number of st.dev. for tracklet cut to keep
		float  dphi        = 0.08,        // dphi window (sigma of tracklet cut)
		float  dtht        = 0.025,       // dtheta .... (if negative, abs will be used with additional cut on |dthetaX|, apart from w.distance
		float  phishift    = 0.0045,      // bending shift
		Bool_t remOvl      = kTRUE,       
		float  injScale    = 1.,//0.7,    // inject injScale*Ncl(Lr1/Lr2) hits
		//
		Bool_t useMC  = kTRUE,           // fill MC info (doRec=kTRUE)
		Float_t etaMin     = -2.,        // min eta range to fill in histos
		Float_t etaMax     =  2.,        // max eta range to fill in histos
		Float_t zMin       = -17,         // process events with Z vertex min
		Float_t zMax       =  17,         //                     max positions
		Float_t scaleMCV0  = 0.8,     // rescale MC V0 to match data
		Float_t ntMin      =   1,         // process events with ESDmult 
		Float_t ntMax      = 20000,       // within this range
		//
		Bool_t checkReconstructables = kFALSE,//kTRUE, // fill histos for reconstructable (needs useMC and doRec) 
		// 
		//---------------------------------------------------------------------------------
		float  phiRot      = 3.14159e+00, // angle for bg. generation with rotation
		float  ovlPhiCut   = 0.005, 
		float  ovlZetaCut  = 0.05,
		TString alirootVer = "VO_ALICE@AliRoot::v4-21-33-AN",
		TString rootVer    = "default",//"VO_ALICE@ROOT::v5-27-06b",
		//
		TString proofCluster="shahoian@alice-caf.cern.ch"
		) 
{ 
  //  
  Bool_t runLocal = kTRUE; // true only for local test mode
  if (runLocal) {
    dataset = "/default/shahoian/tstsim_LHC11d3_146806";
    //dataset = "/default/shahoian/test";
    proofCluster = "";
    alirootVer = "AliRootProofLite";
    nEvents = 500;
  }
  //
  if (!dataset.Contains("sim") && useMC) {
    printf("Running with read data dataset, switching OFF useMC\n");
    useMC = kFALSE;
  }
  //
  printf("Requested: %s %s\n",alirootVer.Data(), rootVer.Data());
  printf("Output expected in %s\n",outFName.Data());
  //
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  //
  TString alirootMode="REC";
  TString extraLibs = "ITSrec:CDB:Geom:"; // not needed in default aliroot mode
  //extraLibs+= "ANALYSIS:ANALYSISalice";
  extraLibs+= "ANALYSIS:OADB:ANALYSISalice:EventMixing";
  TList *list = new TList();
  // sets $ALIROOT_MODE on each worker to let proof to know to run in special mode
  list->Add(new TNamed("ALIROOT_MODE"      , alirootMode.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES", "ITS:include"));
  //  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN","1"));
  //
  //REM: same version of AliRoot on client!!!!! Otherwise error!! 
  TProof::Mgr(proofCluster.Data())->SetROOTVersion(rootVer.Data());
  TProof::Open(proofCluster.Data());//,"workers=10x");
  //  TProof::Open(proofCluster.Data(),"workers=1x");
  if (!gProof) {
    Error("runAAFbgpp.C","Connection to AF failed.");
    return;
  }
  gProof->Exec("TObject *o = gEnv->GetTable()->FindObject(\"Proof.UseMergers\");"
	       "gEnv->GetTable()->Remove(o);", kTRUE);
  //  gProof->SetParameter("PROOF_UseMergers", 0);
  // Lets enable aliroot + extra libs on proof cluster
  if (runLocal) gProof->UploadPackage(alirootVer.Data());
  gProof->EnablePackage(alirootVer.Data(), list);
  //
  gROOT->LoadMacro("MyAnalysisMacroUni.C");

  if (runLocal) {
    Int_t numWorkers = gProof->GetParallel();
    if (numWorkers<1) {printf("No workers\n"); return;}
    gProof->SetParameter("PROOF_PacketizerStrategy", (Int_t)0);
    int frac = (Int_t) 5 / numWorkers;
    if (frac<1) frac = 1;
    gProof->SetParameter("PROOF_PacketAsAFraction", frac);
  }
  MyAnalysisMacroUni(dataset,outFName,noMergeDir,nEvents,useMC,
		     etaMin,etaMax,zMin,zMax,ntMin,ntMax,
		     injScale,scaleDTheta,nStdDev,dphi,dtht,
		     phishift,remOvl,ovlPhiCut,ovlZetaCut,scaleMCV0,checkReconstructables);
  //
}
