/// \file AddTaskRhoDev.C
/// \brief AddTask macro for the AliAnalysisTaskRhoDev class.
///
/// AddTask macro for the AliAnalysisTaskRhoDev class.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date June 16, 2017

AliAnalysisTaskRhoDev* AddTaskRhoDev(
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
 return AliAnalysisTaskRhoDev::AddTaskRhoDev(nTracks, nClusters, nRho, jetradius, acceptance, jetType, rscheme, histo, suffix);
}
