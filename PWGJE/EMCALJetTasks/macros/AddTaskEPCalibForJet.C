AliAnalysisTaskEPCalibForJet* AddTaskEPCalibForJet(TString EPCailbType,
    TString EPCalibJEHandRefFileName,TString EPCalibOrigRefFileName,
    const char *ntracks, const char *nclusters, const char* ncells, const char *suffix)
{
    return AliAnalysisTaskEPCalibForJet::AddTaskEPCalibForJet(
        EPCailbType,
        EPCalibJEHandRefFileName,
        EPCalibOrigRefFileName,
        ntracks, 
        nclusters,
        ncells,
        suffix);
}


