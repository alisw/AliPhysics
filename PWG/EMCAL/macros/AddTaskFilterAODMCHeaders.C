AliEmcalFilterAODMCHeaders* AddTaskFilterAODMCHeaders(
  const char *nParticlesOut = "MCParticlesNotRejected",
  const char *nTracksOut    = "MCTracksNotRejected",
  const char *nClustersOut  = "MCClustersNotRejected",
  Int_t debug               = 0
){
  return AliEmcalFilterAODMCHeaders::AddTaskFilterAODMCHeaders(nParticlesOut, nTracksOut, nClustersOut, debug);
}
