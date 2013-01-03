// $Id$

AliJetResponseMaker* AddTaskJetRespPtHard(const char *ntracks            = "Tracks",
					  const char *nclusters          = "CaloClusters",
					  const char *njets              = "Jets",
					  const char *nmctracks          = "MCParticles",
					  const char *nmcjets            = "MCJets",
					  Double_t    jetradius          = 0.2,
					  Double_t    jetptcut           = 1,
					  Double_t    jetareacut         = 0.557,
					  Double_t    jetBiasTrack       = 5,
					  Double_t    jetBiasClus        = 5,
					  Double_t    maxDistance        = 0.25,
					  UInt_t      type               = AliAnalysisTaskEmcal::kTPC,
					  Int_t       minPtBin           = 1, 
					  Int_t       maxPtBin           = 10,
					  Bool_t      domatch            = kTRUE,
					  Bool_t      biggerMatrix       = kFALSE,
					  const char *taskname           = "AliJetResponseMaker"
					  
)
{  
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C");
  
  AliJetResponseMaker *jetTask = new AliJetResponseMaker[maxPtBin - minPtBin + 1];

  for (Int_t i = minPtBin; i <= maxPtBin; i++) {
    AddTaskJetResponseMaker(ntracks, nclusters, njets, nmctracks, nmcjets,
			    jetradius, jetptcut, jetareacut, jetBiasTrack, 
			    jetBiasClus, maxDistance, type, i, taskname, jetTask + i - minPtBin);
    jetTask[i - minPtBin].SetDoMatching(domatch);
    if (biggerMatrix) 
      jetTask[i - minPtBin].SetHistoBins(1000,0,500);
  }
  
  return jetTask;
}
