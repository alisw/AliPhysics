AliAnalysisTaskScale* AddTaskScale(
  const char *nTracks        = "Tracks",
  const char *nClusters      = "CaloClustersCorr",
  const Double_t trackptcut  = 0.150,
  const Double_t clusptcut   = 0.150,
  const char *taskname       = "Scale",
  const char *sfuncPath      = 0,
  const char *sfuncName      = 0
)
{  
  return AliAnalysisTaskScale::AddTaskScale(nTracks, nClusters, trackptcut, clusptcut, taskname, sfuncPath, sfuncName);
}
