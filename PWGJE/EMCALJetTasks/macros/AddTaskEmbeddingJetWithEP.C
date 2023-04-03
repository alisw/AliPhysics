AliAnalysisTaskEmbeddingJetWithEP* AddTaskEmbeddingJetWithEP(TString EPCailbType,
    TString EPCalibJEHandRefFileName,TString EPCalibOrigRefFileName,
    const char *ntracks, const char *nclusters, const char* ncells, const char *suffix)
{
    return AliAnalysisTaskEmbeddingJetWithEP::AddTaskEmbeddingJetWithEP(
        EPCailbType,
        EPCalibJEHandRefFileName,
        EPCalibOrigRefFileName,
        ntracks, 
        nclusters,
        ncells,
        suffix);
}

