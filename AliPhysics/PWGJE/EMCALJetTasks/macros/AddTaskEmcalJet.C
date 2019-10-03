AliEmcalJetTask* AddTaskEmcalJet(
  const char *nTracks                        = "usedefault",
  const char *nClusters                      = "usedefault",
  const AliJetContainer::EJetAlgo_t jetAlgo  = AliJetContainer::antikt_algorithm,
  const Double_t radius                      = 0.4,
  const AliJetContainer::EJetType_t jetType  = AliJetContainer::kFullJet,
  const Double_t minTrPt                     = 0.15,
  const Double_t minClPt                     = 0.30,
  const Double_t ghostArea                   = 0.005,
  const AliJetContainer::ERecoScheme_t reco  = AliJetContainer::pt_scheme,
  const char *tag                            = "Jet",
  const Double_t minJetPt                    = 0.,
  const Bool_t lockTask                      = kTRUE,
  const Bool_t bFillGhosts                   = kFALSE,
  const char *suffix                         = ""
)
{
  return AliEmcalJetTask::AddTaskEmcalJet(nTracks, nClusters, jetAlgo, radius, jetType, minTrPt, minClPt, ghostArea,
                                          reco, tag, minJetPt, lockTask, bFillGhosts, suffix);
}
