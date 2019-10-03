// $Id$

AliJetResponseMaker* AddTaskJetRespPtHard(  
  const char *ntracks1           = "Tracks",
  const char *nclusters1         = "CaloClusters",
  const char *njets1             = "Jets",
  const char *nrho1              = "Rho",
  Double_t    jetradius1         = 0.2,
  const char *ntracks2           = "MCParticles",
  const char *nclusters2         = "",
  const char *njets2             = "MCJets",
  const char *nrho2              = "",
  Double_t    jetradius2         = 0.2,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.557,
  Double_t    jetBias            = 5,
  Int_t       biasType           = 0,   //  0 = charged, 1 = neutral, 2 = both
  UInt_t      matching           = AliJetResponseMaker::kGeometrical,
  Double_t    maxDistance1       = 0.25,
  Double_t    maxDistance2       = 0.25,
  const char *cutType            = "TPC",
  Int_t       minPtHardBin       = -999,
  Int_t       maxPtHardBin       = -999,
  Int_t       ncent              = 0,
  const char *taskname           = "AliJetResponseMaker",
  Bool_t      biggerMatrix       = kFALSE
)
{
  TCollection *funct = gROOT->GetListOfGlobalFunctions();
  if (!funct->Contains("AddTaskJetResponseMaker"))
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C");
  else
    Printf("Function AddTaskJetResponseMaker already loaded, will not load again...");

  Double_t centRanges[5] = {0,10,30,50,100};
 
  if (ncent == 0) {
    ncent = 1;
    centRanges[0] = -999;
    centRanges[1] = -999;
  }

  if (ncent > 4)
    ncent = 4;

  Int_t ntasks = (maxPtHardBin - minPtHardBin + 1) * ncent;

  if (jetBias > 5)
    ntasks *= 3;
  else if (jetBias > 0)
    ntasks *= 2;

  AliJetResponseMaker *jetTask = new AliJetResponseMaker[ntasks];

  Int_t itask = 0;

  for (Int_t i = minPtHardBin; i <= maxPtHardBin; i++) {
    for (Int_t j = 0; j < ncent; j++) {
      Printf("Adding AliJetResponseMaker n. %d", itask);
      AddTaskJetResponseMaker(ntracks1, nclusters1, njets1, nrho1, jetradius1, ntracks2, nclusters2, njets2, nrho2, jetradius2,
			      jetptcut, jetareacut, 0, biasType,
			      matching, maxDistance1, maxDistance2, cutType, i, centRanges[j], centRanges[j+1], taskname, biggerMatrix, jetTask + itask);
      itask++;

      if (jetBias > 5) {
	Printf("Adding AliJetResponseMaker n. %d", itask);
	AddTaskJetResponseMaker(ntracks1, nclusters1, njets1, nrho1, jetradius1, ntracks2, nclusters2, njets2, nrho2, jetradius2,
				jetptcut, jetareacut, 5, biasType,
				0, 1, 1, cutType, i, centRanges[j], centRanges[j+1], taskname, biggerMatrix, jetTask + itask);
	itask++;
      }

      if (jetBias > 0) {
	Printf("Adding AliJetResponseMaker n. %d", itask);
	AddTaskJetResponseMaker(ntracks1, nclusters1, njets1, nrho1, jetradius1, ntracks2, nclusters2, njets2, nrho2, jetradius2,
				jetptcut, jetareacut, jetBias, biasType,
				0, 1, 1, cutType, i, centRanges[j], centRanges[j+1], taskname, biggerMatrix, jetTask + itask);
	itask++;
      }
    }
  }
  
  return jetTask;
}
