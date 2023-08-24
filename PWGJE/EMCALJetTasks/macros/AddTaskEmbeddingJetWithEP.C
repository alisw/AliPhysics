AliAnalysisTaskEmbeddingJetWithEP* AddTaskEmbeddingJetWithEP(
    const char *ntracks, const char *nclusters, const char* ncells, const char *suffix)
{
    return AliAnalysisTaskEmbeddingJetWithEP::AddTaskEmbeddingJetWithEP(
        ntracks, 
        nclusters,
        ncells,
        suffix);
}

