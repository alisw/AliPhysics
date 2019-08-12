AliAnalysisTaskRho* AddTaskRhoNew (
   const char*    nTracks                        = "usedefault",
   const char*    nClusters                      = "usedefault",
   const char*    nRho                           = "Rho",
   Double_t       jetradius                      = 0.2,
   UInt_t         acceptance                     = AliEmcalJet::kTPCfid,
   AliJetContainer::EJetType_t jetType           = AliJetContainer::kChargedJet,
   const Bool_t   histo                          = kFALSE,
   AliJetContainer::ERecoScheme_t rscheme        = AliJetContainer::pt_scheme,
   const char*    suffix                         = ""
)
{
  return AliAnalysisTaskRho::AddTaskRhoNew(nTracks, nClusters, nRho, jetradius, acceptance, jetType, histo, rscheme, suffix);
}
