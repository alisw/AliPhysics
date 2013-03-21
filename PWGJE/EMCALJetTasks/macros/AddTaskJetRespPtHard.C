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
  UInt_t      matching           = AliJetResponseMaker::kGeometrical,
  Double_t    maxDistance1       = 0.25,
  Double_t    maxDistance2       = 0.25,
  UInt_t      type               = AliAnalysisTaskEmcal::kTPC,
  Int_t       minPtHardBin       = -999,
  Int_t       maxPtHardBin       = -999,
  Int_t       ncent              = 0,
  const char *taskname           = "AliJetResponseMaker",
  Bool_t      biggerMatrix       = kFALSE
)
{  
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C");

  Double_t centRanges[5] = {0,10,30,50,100};
 
  if (ncent == 0) {
    ncent = 1;
    centRanges[0] = -999;
    centRanges[1] = -999;
  }

  if (ncent > 4)
    ncent = 4;

  Int_t ntasks = (maxPtHardBin - minPtHardBin + 1) * ncent;

  if (jetBiasTrack > 5)
    ntasks *= 5;
  else if (jetBiasTrack > 0)
    ntasks *= 5;

  AliJetResponseMaker *jetTask = new AliJetResponseMaker[ntasks];

  Int_t itask = 0;

  for (Int_t i = minPtHardBin; i <= maxPtHardBin; i++) {
    for (Int_t j = 0; j < ncent; j++) {
      Printf("Adding AliJetResponseMaker n. %d", itask);
      AddTaskJetResponseMaker(ntracks1, nclusters1, njets1, nrho1, ntracks2, nclusters2, njets2, nrho2,
			      jetradius, jetptcut, jetareacut, 0, 0,
			      matching, maxDistance1, maxDistance2, type, i, centRanges[j], centRanges[j+1], taskname, biggerMatrix, jetTask + itask);
      itask++;

      if (jetBiasTrack > 5) {
	Printf("Adding AliJetResponseMaker n. %d", itask);
	AddTaskJetResponseMaker(ntracks1, nclusters1, njets1, nrho1, ntracks2, nclusters2, njets2, nrho2,
				jetradius, jetptcut, jetareacut, 5, 1000,
				0, 1, 1, type, i, centRanges[j], centRanges[j+1], taskname, biggerMatrix, jetTask + itask);
	itask++;
      }

      if (jetBiasTrack > 0) {
	Printf("Adding AliJetResponseMaker n. %d", itask);
	AddTaskJetResponseMaker(ntracks1, nclusters1, njets1, nrho1, ntracks2, nclusters2, njets2, nrho2,
				jetradius, jetptcut, jetareacut, jetBiasTrack, 1000,
				0, 1, 1, type, i, centRanges[j], centRanges[j+1], taskname, biggerMatrix, jetTask + itask);
	itask++;
      }
    }
  }
  
  return jetTask;
}
