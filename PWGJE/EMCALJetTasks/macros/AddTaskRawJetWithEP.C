AliAnalysisTaskRawJetWithEP* AddTaskRawJetWithEP(
    const char *ntracks, const char *nclusters, const char* ncells, const char *suffix)
{
    return AliAnalysisTaskRawJetWithEP::AddTaskRawJetWithEP(
        ntracks, 
        nclusters,
        ncells,
        suffix);
}

