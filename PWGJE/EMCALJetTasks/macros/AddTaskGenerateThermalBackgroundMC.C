AliAnalysisTaskGenerateThermalBackgroundMC* AddTaskGenerateThermalBackgroundMC(
  const char *outputCollectionName     = "thermalparticles",
  const Double_t beta                  = 0.3,
  const char *suffix                   = "")
{
  return AliAnalysisTaskGenerateThermalBackgroundMC::AddTaskGenerateThermalBackgroundMC(outputCollectionName,
                                                                                        beta,
                                                                                        suffix);
}
