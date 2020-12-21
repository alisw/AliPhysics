AliAnalysisTaskSDKLResponse* AddTaskSDKLResponse(const char *ntracks = "usedefault",
                                                 const char *njets1 = "Jets1",
                                                 const char *njets2 = "Jets2",
                                                 const char *nrho = "Rho",
                                                 Int_t nCentBins = 1,
                                                 Double_t jetradius = 0.4,
                                                 Double_t jetptcut = 1.0,
                                                 Double_t jetareacut = 0.6,
                                                 const char *type = "EMCAL",
                                                 Int_t backgroption = 0,
                                                 Int_t leadhadtype = 0,
                                                 Double_t fractioneventsfortree = 1.e-6,
                                                 const char *taskname = "AliAnalysisTaskSDKLResponse")
{
  return AliAnalysisTaskSDKLResponse::AddTaskSoftDropResponse(ntracks, njets1, njets2, nrho, nCentBins, jetradius, jetptcut, jetareacut, type, backgroption, leadhadtype, fractioneventsfortree, taskname);
}
