// $Id$

AliJetResponseMaker* AddTaskJetRespPtHard(const char *ntracks            = "Tracks",
					  const char *nclusters          = "CaloClusters",
					  const char *njets              = "Jets",
					  const char *nmcjets            = "MCJets",
					  const char *nmctracks          = "MCParticles",
					  Double_t    jetradius          = 0.4,
					  Double_t    jetptcut           = 1,
					  Double_t    jetareacut         = 0.8,
					  Double_t    ptcut              = 0.15,
					  Double_t    jetBiasTrack       = 10,
					  Double_t    jetBiasClus        = 10,
					  Double_t    maxDistance        = 0.25,
					  UInt_t      type               = AliAnalysisTaskEmcal::kTPC,
					  Int_t minPtBin                 = 1, 
					  Int_t maxPtBin                 = 11,
					  const char *taskname           = "AliJetResponseMaker"
)
{  
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/AddTaskJetResponseMaker.C");
  
  AliJetResponseMaker *jetTask = 0;

  for (Int_t i = minPtBin; i <= maxPtBin; i++) {
    jetTask = AddTaskJetResponseMaker(ntracks, nclusters, njets, nmcjets, nmctracks, 
				      jetradius, jetptcut, jetareacut, ptcut, jetBiasTrack, 
				      jetBiasClus, maxDistance, type, i, taskname);
  }
  
  return jetTask;
}
