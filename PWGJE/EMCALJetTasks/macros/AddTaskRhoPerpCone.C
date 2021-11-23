/// \file AddTaskRhoPerpCone.C
/// \brief AddTask macro for the AliAnalysisTaskRhoPerpCone class.
///
/// AddTask macro for the AliAnalysisTaskRhoPerpCone class.
///

AliAnalysisTaskRhoPerpCone* AddTaskRhoPerpCone(
    TString        nTracks                        = "usedefault",
    Double_t       trackPtCut                     = 0.15,
    TString        nClusters                      = "usedefault",
    Double_t       clusECut                       = 0.30,
    TString        nRho                           = "Rho",
    Double_t       jetradius                      = 0.2,
    UInt_t         acceptance                     = AliEmcalJet::kTPCfid,
    AliJetContainer::EJetType_t jetType           = AliJetContainer::kChargedJet,
    AliJetContainer::ERecoScheme_t rscheme        = AliJetContainer::pt_scheme,
    Bool_t         histo                          = kTRUE,
    TString        suffix                         = ""
)
{  
 return AliAnalysisTaskRhoPerpCone::AddTaskRhoPerpCone(nTracks, trackPtCut, nClusters, clusECut, nRho, jetradius, acceptance, jetType, rscheme, histo, suffix);
}
