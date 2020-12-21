/// \file AddTaskRhoBaseDev.C
/// \brief AddTask macro for the AliAnalysisTaskRhoBaseDev class.
///
/// AddTask macro for the AliAnalysisTaskRhoBaseDev class.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date June 16, 2017

AliAnalysisTaskRhoBaseDev* AddTaskRhoBaseDev(
    TString        nTracks                        = "usedefault",
    TString        nClusters                      = "usedefault",
    TString        nRho                           = "Rho",
    Double_t       jetradius                      = 0.2,
    UInt_t         acceptance                     = AliEmcalJet::kTPCfid,
    AliJetContainer::EJetType_t jetType           = AliJetContainer::kChargedJet,
    AliJetContainer::ERecoScheme_t rscheme        = AliJetContainer::pt_scheme,
    Bool_t         histo                          = kTRUE,
    TString        suffix                         = ""
)
{  
 return AliAnalysisTaskRhoBaseDev::AddTaskRhoBaseDev(nTracks, nClusters, nRho, jetradius, acceptance, jetType, rscheme, histo, suffix);
}
