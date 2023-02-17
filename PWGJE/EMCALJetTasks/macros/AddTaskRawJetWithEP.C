AliAnalysisTaskRawJetWithEP* AddTaskRawJetWithEP(TString EPCailbType,
    TString EPCalibJEHandRefFileName,TString EPCalibOrigRefFileName,
    const char *ntracks, const char *nclusters, const char* ncells, const char *suffix)
{
    return AliAnalysisTaskRawJetWithEP::AddTaskRawJetWithEP(
        EPCailbType,
        EPCalibJEHandRefFileName,
        EPCalibOrigRefFileName,
        ntracks, 
        nclusters,
        ncells, 
        suffix);
}

