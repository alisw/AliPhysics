AliAnalysisTaskSDKLResponse* AddTaskSDKLResponse(const char *ntracks, const char *njets1, const char *njets2, const char *nrho, Int_t nCentBins,
                                                                       Double_t jetradius, Double_t jetptcut, Double_t jetareacut, const char *type, Int_t backgroption,
                                                                       Int_t leadhadtype, Double_t fractioneventsfortree, const char *taskname)
{
  return AliAnalysisTaskSDKLResponse::AddTaskSoftDropResponse(ntracks, njets1, njets2, nrho, nCentBins, jetradius, jetptcut, jetareacut, type, backgroption, leadhadtype, fractioneventsfortree, taskname);
}