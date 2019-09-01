AliAnalysisTaskSDKL* AddTaskSDKL(const char *ntracks = "usedefault",
                                 const char *njets = "Jets",
                                 const char *nrho = "Rho",
                                 Int_t nCentBins = 1,
                                 Double_t jetradius = 0.4,
                                 Double_t jetptcut = 1.0,
                                 Double_t jetareacut = 0.6,
                                 const char *type = "EMCAL",
                                 Int_t backgroption = 0,
                                 Int_t leadhadtype = 0,
                                 const char *taskname = "AliAnalysisTaskSDKL")
{

   return AliAnalysisTaskSDKL::AddTaskSoftDrop(ntracks, njets, nrho, nCentBins, jetradius, jetptcut, jetareacut, type, backgroption, leadhadtype, taskname);

}
