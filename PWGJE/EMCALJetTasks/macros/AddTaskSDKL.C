AliAnalysisTaskSDKL* AddTaskSDKL(const char *ntracks, const char *njets, const char *nrho, Int_t nCentBins, Double_t jetradius, 
                                 Double_t jetptcut, Double_t jetareacut, const char *type, Int_t backgroption, 
                                 Int_t leadhadtype, const char *taskname)
{

   return AliAnalysisTaskSDKL::AddTaskSoftDrop(ntracks, njets, nrho, nCentBins, jetradius, jetptcut, jetareacut, type, backgroption, leadhadtype, taskname);

}