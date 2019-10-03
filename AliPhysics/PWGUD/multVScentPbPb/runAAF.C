//
void runAAF(TString dataset="/alice/sim/LHC10f8c_130844",
	    TString outFName = "trbg.root",
	    Bool_t doRec  = kTRUE,           // fill data histos from new reco
	    Bool_t doInj  = kTRUE,           // create Inj. bg
	    Bool_t doRot  = kTRUE,           // create Rot. bg
	    Bool_t doMix  = kFALSE,//kTRUE,  // create Mix. bg
	    Bool_t useMC  = kTRUE,           // fill MC info (doRec=kTRUE)
	    Bool_t checkReconstructables = kFALSE,//kTRUE, // fill histos for reconstructable (needs useMC and doRec) 
	    // 
	    Float_t etaCut     = 3.0,        // max |eta| range to fill in histos
	    //
	    // specific parameters for reconstruction
	    //----------------------- Zv selection parameters important for mixing, to be tuned
	    Float_t zMin       = -20,        // process events with Z vertex min
	    Float_t zMax       =  20,        //                     max positions
	    Float_t zMixBinSz  =  20, //0.1,  // Zv. bin for mixing
	    //---------------------------------------------------------------------------------
	    //
	    //----------------------- Ntracklets selection parameters important for mixing, to be tuned
	    Float_t ntMin      =   1,         // process events with ESDmult 
	    Float_t ntMax      = 20000,       // within this range
	    Float_t ntMixBinSz = 20000,       // ESDMult bin size for mixing
	    //---------------------------------------------------------------------------------
	    float  phiRot      = 3.14159e+00, // angle for bg. generation with rotation
	    float  injScale    = 1.,//0.7,    // inject injScale*Ncl(Lr1/Lr2) hits
	    Bool_t scaleDTheta = kTRUE,       // scale dTheta by 1/sin^2(theta) in trackleting
	    float  nStdDev     = 25.,         // number of st.dev. for tracklet cut to keep
	    float  dphi        = 0.06,        // dphi window (sigma of tracklet cut)
	    float  dtht        = 0.025,       // dtheta .... (if negative, abs will be used with additional cut on |dthetaX|, apart from w.distance
	    float  phishift    = 0.0045,      // bending shift
	    Bool_t remOvl      = kTRUE,       
	    float  ovlPhiCut   = 0.005, 
	    float  ovlZetaCut  = 0.05,
            Int_t  nEvents     = 50000000,
	    Int_t  nEventsSkip = 0,
	    TString alirootVer = "VO_ALICE@AliRoot::v4-21-05-AN",
	    TString rootVer    = "default",//"VO_ALICE@ROOT::v5-27-06b",
	    //
	    TString proofCluster="shahoian@alice-caf.cern.ch"
	    ) 
{ 
  //  
  Bool_t runLocal = kFALSE; // true only for local test mode
  if (runLocal) {
    dataset = "/default/shahoian/test_pp";
    //dataset = "/default/shahoian/test";
    proofCluster = "";
    alirootVer = "AliRootProofLite";
    nEvents = 500;
  }
  //
  if (dataset.Contains("alice/data") && useMC) {
    printf("Running with read data dataset, switching OFF useMC\n");
    useMC = kFALSE;
  }
  //
  printf("Start Analysis for %s, max %d Events skipping %d\n"
	 "Event Cuts: |eta|<%.1f, %.2f<Zv<%.2f (Mix.Bin:%.2f), %.0f<Mult<%0.f (Mix.Bin:%.0f)\n"
	 "Reco:%d Inj:%d Rot:%d Mix:%d | MCinfo:%d CheckReconstructables:%d\n"
	 "PhiRot:%.4f InjScale:%.2f ScaleDTheta:%d NStdDev:%.1f DPhi:%.4f DTheta:%.4f PhiShift:%.4f\n"
	 "RemoveOverlaps:%d PhiOvl:%.4f ZEtaOvl:%.4f\n",
	 dataset.Data(),nEvents,nEventsSkip,
	 etaCut,zMin,zMax,zMixBinSz,ntMin,ntMax,ntMixBinSz,
	 doRec,doInj,doRot,doMix,useMC,checkReconstructables,
	 phiRot,injScale,scaleDTheta,nStdDev,dphi,dtht,phishift,
	 remOvl,ovlPhiCut,ovlZetaCut);
  //
  printf("Requested: %s %s\n",alirootVer.Data(), rootVer.Data());
  printf("Output expected in %s\n",outFName.Data());
  //
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  //
  TString alirootMode="REC";
  TString extraLibs = "ITSrec:CDB:Geom:"; // not needed in default aliroot mode
  //extraLibs+= "ANALYSIS:ANALYSISalice";
  extraLibs+= "ANALYSIS:ANALYSISalice:EventMixing";
  TList *list = new TList();
  // sets $ALIROOT_MODE on each worker to let proof to know to run in special mode
  list->Add(new TNamed("ALIROOT_MODE"      , alirootMode.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN","1"));
  //
  //REM: same version of AliRoot on client!!!!! Otherwise error!! 
  TProof::Mgr(proofCluster.Data())->SetROOTVersion(rootVer.Data());
  TProof::Open(proofCluster.Data());//,"workers=10x");
  //  TProof::Open(proofCluster.Data(),"workers=1x");
  if (!gProof) {
    Error("runAAF.C","Connection to AF failed.");
    return;
  }
  gProof->Exec("TObject *o = gEnv->GetTable()->FindObject(\"Proof.UseMergers\");"
	       "gEnv->GetTable()->Remove(o);", kTRUE);
  //  gProof->SetParameter("PROOF_UseMergers", 0);
  // Lets enable aliroot + extra libs on proof cluster
  if (runLocal) gProof->UploadPackage(alirootVer.Data());
  gProof->EnablePackage(alirootVer.Data(), list);
  //
  gROOT->LoadMacro("MyAnalysisMacro.C");

  if (runLocal) {
    Int_t numWorkers = gProof->GetParallel();
    if (numWorkers<1) {printf("No workers\n"); return;}
    gProof->SetParameter("PROOF_PacketizerStrategy", (Int_t)0);
    int frac = (Int_t) 5 / numWorkers;
    if (frac<1) frac = 1;
    gProof->SetParameter("PROOF_PacketAsAFraction", frac);
  }
  MyAnalysisMacro(dataset,outFName,doRec,doInj,doRot,doMix,useMC,checkReconstructables,
		  etaCut,zMin,zMax,zMixBinSz,ntMin,ntMax,ntMixBinSz,
		  phiRot,injScale,scaleDTheta,nStdDev,dphi,dtht,
		  phishift,remOvl,ovlPhiCut,ovlZetaCut,nEvents,nEventsSkip);
  //
}
