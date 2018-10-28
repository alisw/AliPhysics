// $Id$

AliAnalysisTaskRhoSparse* AddTaskRhoSparse(
					   const char    *nTracks     = "usedefault",
					   const char    *nClusters   = "usedefault",
					   const char    *nRho        = "Rho",
					   Double_t       jetradius   = 0.2,
					   UInt_t         acceptance  = AliEmcalJet::kTPCfid,
					   AliJetContainer::EJetType_t jetType    = AliJetContainer::kChargedJet,
					   AliJetContainer::ERecoScheme_t rscheme = AliJetContainer::pt_scheme,
					   const Bool_t   histo       = kFALSE,
					   const char    *nJetsSig    = "",
					   const char    *cutType     = "TPC",
					   Double_t       jetptcut    = 0.0,
					   Double_t       jetareacut  = 0.01,
					   Double_t       emcareacut  = 0,
					   const char    *suffix      = ""
					   )
{  
	return AliAnalysisTaskRhoSparse::AddTaskRhoSparse(nTracks,
													 nClusters,
													 nRho,
													 jetradius,
													 acceptance,
													 jetType,
													 rscheme,
													 histo,
													 nJetsSig,
													 cutType,
													 jetptcut,
													 jetareacut,
													 emcareacut,
													 suffix);
}
