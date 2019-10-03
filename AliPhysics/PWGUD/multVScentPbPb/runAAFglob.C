//
void runAAFglob(TString dataset="/alice/data/LHC10h_000137161_p1_plusplusplus",
		TString outFName = "globs.root",
		Int_t  nEvents     = 50000000,
		Int_t  nEventsSkip = 0,
		TString alirootVer = "VO_ALICE@AliRoot::v4-21-04-AN",
		TString rootVer    = "VO_ALICE@ROOT::v5-27-06b",
		TString proofCluster="shahoian@alice-caf.cern.ch"
		//TString proofCluster="shahoian@skaf.saske.sk"
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
  printf("Requested: %s %s\n",alirootVer.Data(), rootVer.Data());
  printf("Output expected in %s\n",outFName.Data());
  //
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  //
  TString alirootMode="REC";
  TString extraLibs;
  extraLibs+= "ANALYSIS:ANALYSISalice";
  TList *list = new TList();
  // sets $ALIROOT_MODE on each worker to let proof to know to run in special mode
  list->Add(new TNamed("ALIROOT_MODE"      , alirootMode.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
  //  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN","1"));
  //
  //REM: same version of AliRoot on client!!!!! Otherwise error!! 
  TProof::Mgr(proofCluster.Data())->SetROOTVersion(rootVer.Data());
  TProof::Open(proofCluster.Data());//,"workers=10x");
  //  TProof::Open(proofCluster.Data(),"workers=1x");
  if (!gProof) {
    Error("runAAFglob.C","Connection to AF failed.");
    return;
  }
  gProof->Exec("TObject *o = gEnv->GetTable()->FindObject(\"Proof.UseMergers\");"
	       "gEnv->GetTable()->Remove(o);", kTRUE);
  //  gProof->SetParameter("PROOF_UseMergers", 0);
  // Lets enable aliroot + extra libs on proof cluster
  if (runLocal) gProof->UploadPackage(alirootVer.Data());
  gProof->EnablePackage(alirootVer.Data(), list);
  //
  gROOT->LoadMacro("AnalysisMacroGlob.C");
  //
  if (runLocal) {
    Int_t numWorkers = gProof->GetParallel();
    if (numWorkers<1) {printf("No workers\n"); return;}
    gProof->SetParameter("PROOF_PacketizerStrategy", (Int_t)0);
    int frac = (Int_t) 5 / numWorkers;
    if (frac<1) frac = 1;
    gProof->SetParameter("PROOF_PacketAsAFraction", frac);
  }
  AnalysisMacroGlob(dataset,outFName,nEvents,nEventsSkip);
  //
}
