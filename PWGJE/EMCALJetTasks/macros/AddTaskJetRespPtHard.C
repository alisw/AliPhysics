// $Id$

AliJetResponseMaker* AddTaskJetRespPtHard(  
  const char *ntracks1           = "Tracks",
  const char *nclusters1         = "CaloClusters",
  const char *njets1             = "Jets",
  const char *nrho1              = "Rho",
  const char *ntracks2           = "MCParticles",
  const char *nclusters2         = "",
  const char *njets2             = "MCJets",
  const char *nrho2              = "",
  Double_t    jetradius          = 0.2,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.557,
  Double_t    jetBiasTrack       = 5,
  Double_t    jetBiasClus        = 5,
  UInt_t      matching           = AliJetResponseMaker::kGeometrical,
  Double_t    maxDistance1       = 0.25,
  Double_t    maxDistance2       = 0.25,
  UInt_t      type               = AliAnalysisTaskEmcal::kTPC,
  Int_t       maxPtHardBin       = -999,
  Int_t       minPtHardBin       = -999,
  const char *taskname           = "AliJetResponseMaker",
  Bool_t      biggerMatrix       = kFALSE
)
{  
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C");
  
  AliJetResponseMaker *jetTask = new AliJetResponseMaker[maxPtHardBin - minPtHardBin + 1];

  for (Int_t i = minPtHardBin; i <= maxPtHardBin; i++) {
    AddTaskJetResponseMaker(ntracks1, nclusters1, njets1, nrho1, ntracks2, nclusters2, njets2, nrho2,
			    jetradius, jetptcut, jetareacut, jetBiasTrack, jetBiasClus, 
			    matching, maxDistance1, maxDistance2, type, i, taskname, biggerMatrix, jetTask + i - minPtHardBin);
  }
  
  return jetTask;
}
