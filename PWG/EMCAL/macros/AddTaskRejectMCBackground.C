AliEmcalRejectMCBackground* AddTaskRejectMCBackground(
  const char *nParticlesOut = "MCParticlesNotRejected",
  const char *nTracksOut    = "MCTracksNotRejected",
  const char *nClustersOut  = "MCClustersNotRejected",
  Int_t signalRejection     = 2,
  Int_t debug               = 0
){
  return AliEmcalRejectMCBackground::AddTaskRejectMCBackground(nParticlesOut, nTracksOut, nClustersOut, signalRejection, debug);
}
