//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis
//  - MC output
//  implementation file
//
//*-- Author: Marcelo G. Munhoz (USP)
	      //_________________________________________________________________________

#include "AliAnalysisEmEtMonteCarlo.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliTrackReference.h"
#include "AliESDEvent.h"
#include "TH2F.h"
#include "TParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "TList.h"
#include "AliESDCaloCluster.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
#include "AliEMCALTrack.h"
#include "AliESDtrackCuts.h"
#include "AliEMCALGeometry.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"
#include "TGeoManager.h"

  using namespace std;

ClassImp(AliAnalysisEmEtMonteCarlo);


// ctor
AliAnalysisEmEtMonteCarlo::AliAnalysisEmEtMonteCarlo():AliAnalysisEtMonteCarlo()
						      ,fPrimtotET(0), fPrimAcctotET(0), fPrimRectotET(0), fPrimRectotETDep(0)
						      ,fElectrontotET(0), fElectronAcctotET(0), fElectronRectotET(0)
						      ,fConvElectrontotET(0), fConvElectronAcctotET(0), fConvElectronRectotET(0), fScatElectrontotET(0), fScatElectronAcctotET(0), fScatElectronRectotET(0)
						      ,fTotElectrontotET(0), fTotElectronAcctotET(0), fTotElectronRectotET(0)

						      ,fGammatotET(0), fGammaAcctotET(0), fGammaRectotET(0)
						      ,fAnnihGammatotET(0), fAnnihGammaAcctotET(0), fAnnihGammaRectotET(0), fScatGammatotET(0), fScatGammaAcctotET(0), fScatGammaRectotET(0)
						      ,fTotGammatotET(0), fTotGammaAcctotET(0), fTotGammaRectotET(0)
						      ,fConvGammatotET(0),fNonConvGammatotET(0),fConvGammaAcctotET(0),fNonConvGammaAcctotET(0), fNPPPi0GammatotET(0), fNPPPi0GammaRectotET(0)

						      ,fTotEMtotET(0), fTotEMAcctotET(0), fTotEMRectotET(0)

						      ,fNPPElectrontotET(0), fNPPElectronRectotET(0), fNPPGammatotET(0), fNPPGammaRectotET(0)
						      ,fTotNPPEMtotET(0), fTotNPPEMRectotET(0)

						      ,fMuontotET(0), fPiontotET(0), fKaontotET(0), fProtontotET(0)
						      ,fMuonAcctotET(0), fPionAcctotET(0), fKaonAcctotET(0), fProtonAcctotET(0)
						      ,fMuonRectotET(0), fMuonRectotETDep(0), fPionRectotET(0), fPionRectotETDep(0), fKaonRectotET(0), fKaonRectotETDep(0), fProtonRectotET(0), fProtonRectotETDep(0)
						      ,fMuonMatchtotET(0), fMuonMatchtotETDep(0), fPionMatchtotET(0), fPionMatchtotETDep(0), fKaonMatchtotET(0), fKaonMatchtotETDep(0), fProtonMatchtotET(0), fProtonMatchtotETDep(0)
						      ,fTotChargedtotET(0), fTotChargedAcctotET(0), fTotChargedRectotET(0), fTotChargedRectotETDep(0), fTotChargedMatchtotET(0), fTotChargedMatchtotETDep(0)

						      ,fNeutrontotET(0), fNeutronAcctotET(0), fNeutronRectotET(0), fNeutronRectotETDep(0)
						      ,fK0totET(0), fK0RectotET(0), fK0RectotETDep(0), fLambdatotET(0), fLambdaRectotET(0), fLambdaRectotETDep(0)
						      ,fTotNeutraltotET(0), fTotNeutralRectotET(0), fTotNeutralRectotETDep(0)

						      ,fTotaltotET(0), fTotalAcctotET(0), fTotalRectotET(0), fTotalRectotETDep(0)

 						      ,fGeoUt(0)

						      ,fHistPrimEtaEET(0) 
						      ,fHistPrimEtaPtET(0) 
						      ,fHistPrimEtaET(0) 
						      ,fHistPrimtotET(0) 

						      ,fHistPrimAccEtaEET(0) 
						      ,fHistPrimAccEtaPtET(0) 
						      ,fHistPrimAccEtaET(0) 
						      ,fHistPrimAcctotET(0) 

						      ,fHistPrimRecEtaEET(0) 
						      ,fHistPrimRecEtaPtET(0) 
						      ,fHistPrimRecEtaET(0) 
						      ,fHistPrimRectotET(0) 

						      ,fHistPrimRecEtaEDepETDep(0) 
						      ,fHistPrimRecEtaPtETDep(0) 
						      ,fHistPrimRecEtaETDep(0) 
						      ,fHistPrimRectotETDep(0) 

						      ,fHistElectronEtaEET(0) 
						      ,fHistElectronEtaPtET(0) 
						      ,fHistElectronEtaET(0) 
						      ,fHistElectronEtaE(0) 
						      ,fHistElectronEtaPt(0) 
						      ,fHistElectrontotET(0) 

						      ,fHistConvElectronEtaEET(0)  
						      ,fHistConvElectronEtaPtET(0)  
						      ,fHistConvElectronEtaET(0)  
						      ,fHistConvElectronEtaE(0)  
						      ,fHistConvElectronEtaPt(0)  
						      ,fHistConvElectrontotET(0)  

						      ,fHistScatElectronEtaEET(0)  
						      ,fHistScatElectronEtaPtET(0)  
						      ,fHistScatElectronEtaET(0)  
						      ,fHistScatElectronEtaE(0)  
						      ,fHistScatElectronEtaPt(0)  
						      ,fHistScatElectrontotET(0)  

						      ,fHistTotElectrontotET(0)

						      ,fHistGammaEtaEET(0)  
						      ,fHistGammaEtaPtET(0)  
						      ,fHistGammaEtaET(0)  
						      ,fHistGammaEtaE(0)  
						      ,fHistGammaEtaPt(0)  
						      ,fHistGammatotET(0)  

						      ,fHistAnnihGammaEtaEET(0)  
						      ,fHistAnnihGammaEtaPtET(0)  
						      ,fHistAnnihGammaEtaET(0)  
						      ,fHistAnnihGammaEtaE(0)  
						      ,fHistAnnihGammaEtaPt(0)  
						      ,fHistAnnihGammatotET(0)  

						      ,fHistScatGammaEtaEET(0)  
						      ,fHistScatGammaEtaPtET(0)  
						      ,fHistScatGammaEtaET(0)  
						      ,fHistScatGammaEtaE(0)  
						      ,fHistScatGammaEtaPt(0)  
						      ,fHistScatGammatotET(0)  

						      ,fHistConvGammaEtaEET(0)  
						      ,fHistConvGammaEtaPtET(0)  
						      ,fHistConvGammaEtaET(0)  
						      ,fHistConvGammaEtaE(0)  
						      ,fHistConvGammaEtaPt(0)  
						      ,fHistConvGammatotET(0)  

						      ,fHistNonConvGammaEtaEET(0)  
						      ,fHistNonConvGammaEtaPtET(0)  
						      ,fHistNonConvGammaEtaET(0)  
						      ,fHistNonConvGammaEtaE(0)  
						      ,fHistNonConvGammaEtaPt(0)  
						      ,fHistNonConvGammatotET(0)  

						      ,fHistTotGammatotET(0)

						      ,fHistTotEMtotET(0)

						      ,fHistNPPElectronEtaEET(0) 
						      ,fHistNPPElectronEtaPtET(0) 
						      ,fHistNPPElectronEtaET(0) 
						      ,fHistNPPElectronEtaE(0) 
						      ,fHistNPPElectronEtaPt(0) 
						      ,fHistNPPElectrontotET(0) 

						      ,fHistNPPGammaEtaEET(0) 
						      ,fHistNPPGammaEtaPtET(0) 
						      ,fHistNPPGammaEtaET(0) 
						      ,fHistNPPGammaEtaE(0) 
						      ,fHistNPPGammaEtaPt(0) 
						      ,fHistNPPGammatotET(0) 

						      ,fHistTotNPPEMtotET(0)

						      ,fHistNPPPi0GammaEtaEET(0) 
						      ,fHistNPPPi0GammaEtaPtET(0) 
						      ,fHistNPPPi0GammaEtaET(0) 
						      ,fHistNPPPi0GammaEtaE(0) 
						      ,fHistNPPPi0GammaEtaPt(0) 
						      ,fHistNPPPi0GammatotET(0) 

						      ,fHistElectronAccEtaEET(0) 
						      ,fHistElectronAccEtaPtET(0) 
						      ,fHistElectronAccEtaET(0) 
						      ,fHistElectronAccEtaE(0) 
						      ,fHistElectronAccEtaPt(0) 
						      ,fHistElectronAcctotET(0) 

						      ,fHistConvElectronAccEtaEET(0)  
						      ,fHistConvElectronAccEtaPtET(0)  
						      ,fHistConvElectronAccEtaET(0)  
						      ,fHistConvElectronAccEtaE(0)  
						      ,fHistConvElectronAccEtaPt(0)  
						      ,fHistConvElectronAcctotET(0)  

						      ,fHistScatElectronAccEtaEET(0)  
						      ,fHistScatElectronAccEtaPtET(0)  
						      ,fHistScatElectronAccEtaET(0)  
						      ,fHistScatElectronAccEtaE(0)  
						      ,fHistScatElectronAccEtaPt(0)  
						      ,fHistScatElectronAcctotET(0)  

						      ,fHistTotElectronAcctotET(0)

						      ,fHistGammaAccEtaEET(0)  
						      ,fHistGammaAccEtaPtET(0)  
						      ,fHistGammaAccEtaET(0)  
						      ,fHistGammaAccEtaE(0)  
						      ,fHistGammaAccEtaPt(0)  
						      ,fHistGammaAcctotET(0)  

						      ,fHistAnnihGammaAccEtaEET(0)  
						      ,fHistAnnihGammaAccEtaPtET(0)  
						      ,fHistAnnihGammaAccEtaET(0)  
						      ,fHistAnnihGammaAccEtaE(0)  
						      ,fHistAnnihGammaAccEtaPt(0)  
						      ,fHistAnnihGammaAcctotET(0)  

						      ,fHistScatGammaAccEtaEET(0)  
						      ,fHistScatGammaAccEtaPtET(0)  
						      ,fHistScatGammaAccEtaET(0)  
						      ,fHistScatGammaAccEtaE(0)  
						      ,fHistScatGammaAccEtaPt(0)  
						      ,fHistScatGammaAcctotET(0)  

						      ,fHistConvGammaAccEtaEET(0)  
						      ,fHistConvGammaAccEtaPtET(0)  
						      ,fHistConvGammaAccEtaET(0)  
						      ,fHistConvGammaAccEtaE(0)  
						      ,fHistConvGammaAccEtaPt(0)  
						      ,fHistConvGammaAcctotET(0)  

						      ,fHistNonConvGammaAccEtaEET(0)  
						      ,fHistNonConvGammaAccEtaPtET(0)  
						      ,fHistNonConvGammaAccEtaET(0)  
						      ,fHistNonConvGammaAccEtaE(0)  
						      ,fHistNonConvGammaAccEtaPt(0)  
						      ,fHistNonConvGammaAcctotET(0)  

						      ,fHistTotGammaAcctotET(0)

						      ,fHistTotEMAcctotET(0)

						      ,fHistNPPElectronAccEtaEET(0) 
						      ,fHistNPPElectronAccEtaPtET(0) 
						      ,fHistNPPElectronAccEtaE(0) 
						      ,fHistNPPElectronAccEtaPt(0) 

						      ,fHistNPPGammaAccEtaEET(0) 
						      ,fHistNPPGammaAccEtaPtET(0) 
						      ,fHistNPPGammaAccEtaE(0) 
						      ,fHistNPPGammaAccEtaPt(0) 

						      ,fHistElectronRecEtaEET(0) 
						      ,fHistElectronRecEtaPtET(0) 
						      ,fHistElectronRecEtaET(0) 
						      ,fHistElectronRecEtaE(0) 
						      ,fHistElectronRecEtaPt(0) 
						      ,fHistElectronRectotET(0) 

						      ,fHistConvElectronRecEtaEET(0)  
						      ,fHistConvElectronRecEtaPtET(0)  
						      ,fHistConvElectronRecEtaET(0)  
						      ,fHistConvElectronRecEtaE(0)  
						      ,fHistConvElectronRecEtaPt(0)  
						      ,fHistConvElectronRectotET(0)  

						      ,fHistScatElectronRecEtaEET(0)  
						      ,fHistScatElectronRecEtaPtET(0)  
						      ,fHistScatElectronRecEtaET(0)  
						      ,fHistScatElectronRecEtaE(0)  
						      ,fHistScatElectronRecEtaPt(0)  
						      ,fHistScatElectronRectotET(0)  

						      ,fHistTotElectronRectotET(0)

						      ,fHistGammaRecEtaEET(0)  
						      ,fHistGammaRecEtaPtET(0)  
						      ,fHistGammaRecEtaET(0)  
						      ,fHistGammaRecEtaE(0)  
						      ,fHistGammaRecEtaPt(0)  
						      ,fHistGammaRectotET(0)  

						      ,fHistAnnihGammaRecEtaEET(0)  
						      ,fHistAnnihGammaRecEtaPtET(0)  
						      ,fHistAnnihGammaRecEtaET(0)  
						      ,fHistAnnihGammaRecEtaE(0)  
						      ,fHistAnnihGammaRecEtaPt(0)  
						      ,fHistAnnihGammaRectotET(0)  

						      ,fHistScatGammaRecEtaEET(0)  
						      ,fHistScatGammaRecEtaPtET(0)  
						      ,fHistScatGammaRecEtaET(0)  
						      ,fHistScatGammaRecEtaE(0)  
						      ,fHistScatGammaRecEtaPt(0)  
						      ,fHistScatGammaRectotET(0)  

						      ,fHistTotGammaRectotET(0)

						      ,fHistTotEMRectotET(0)

						      ,fHistNPPElectronRecEtaEET(0) 
						      ,fHistNPPElectronRecEtaPtET(0) 
						      ,fHistNPPElectronRecEtaET(0) 
						      ,fHistNPPElectronRecEtaE(0) 
						      ,fHistNPPElectronRecEtaPt(0) 
						      ,fHistNPPElectronRectotET(0) 

						      ,fHistNPPGammaRecEtaEET(0) 
						      ,fHistNPPGammaRecEtaPtET(0) 
						      ,fHistNPPGammaRecEtaET(0) 
						      ,fHistNPPGammaRecEtaE(0) 
						      ,fHistNPPGammaRecEtaPt(0) 
						      ,fHistNPPGammaRectotET(0) 

						      ,fHistTotNPPEMRectotET(0)

						      ,fHistNPPPi0GammaRecEtaEET(0) 
						      ,fHistNPPPi0GammaRecEtaPtET(0) 
						      ,fHistNPPPi0GammaRecEtaET(0) 
						      ,fHistNPPPi0GammaRecEtaE(0) 
						      ,fHistNPPPi0GammaRecEtaPt(0) 
						      ,fHistNPPPi0GammaRectotET(0) 

						      ,fHistMuonEtaEET(0) 
						      ,fHistMuonAccEtaEET(0) 
						      ,fHistMuonRecEtaEET(0) 
						      ,fHistMuonMatchEtaEET(0) 

						      ,fHistMuonEtaPtET(0) 
						      ,fHistMuonAccEtaPtET(0) 
						      ,fHistMuonRecEtaPtET(0) 
						      ,fHistMuonMatchEtaPtET(0) 

						      ,fHistMuonEtaET(0) 
						      ,fHistMuonAccEtaET(0) 
						      ,fHistMuonRecEtaET(0) 
						      ,fHistMuonMatchEtaET(0) 

						      ,fHistMuonEtaE(0) 
						      ,fHistMuonAccEtaE(0) 
						      ,fHistMuonRecEtaE(0) 
						      ,fHistMuonMatchEtaE(0) 

						      ,fHistMuonEtaPt(0) 
						      ,fHistMuonAccEtaPt(0) 
						      ,fHistMuonRecEtaPt(0) 
						      ,fHistMuonMatchEtaPt(0) 

						      ,fHistMuontotET(0) 
						      ,fHistMuonAcctotET(0) 
						      ,fHistMuonRectotET(0) 
						      ,fHistMuonMatchtotET(0) 

						      ,fHistMuonRectotETDep(0) 
						      ,fHistMuonMatchtotETDep(0) 

						      ,fHistMuonRecEtaEDepETDep(0) 
						      ,fHistMuonMatchEtaEDepETDep(0) 

						      ,fHistMuonRecEtaPtETDep(0) 
						      ,fHistMuonMatchEtaPtETDep(0) 

						      ,fHistMuonRecEtaETDep(0) 
						      ,fHistMuonMatchEtaETDep(0) 

						      ,fHistMuonRecResEET(0) 
						      ,fHistMuonRecResPtET(0) 
						      ,fHistMuonRecResE(0) 
						      ,fHistMuonRecResPt(0) 

						      ,fHistMuonRecResEDepETDep(0) 
						      ,fHistMuonRecResPtETDep(0) 

						      ,fHistPionEtaEET(0) 
						      ,fHistPionAccEtaEET(0) 
						      ,fHistPionRecEtaEET(0) 
						      ,fHistPionMatchEtaEET(0) 

						      ,fHistPionEtaPtET(0) 
						      ,fHistPionAccEtaPtET(0) 
						      ,fHistPionRecEtaPtET(0) 
						      ,fHistPionMatchEtaPtET(0) 

						      ,fHistPionEtaET(0) 
						      ,fHistPionAccEtaET(0) 
						      ,fHistPionRecEtaET(0) 
						      ,fHistPionMatchEtaET(0) 

						      ,fHistPionEtaE(0) 
						      ,fHistPionAccEtaE(0) 
						      ,fHistPionRecEtaE(0) 
						      ,fHistPionMatchEtaE(0) 

						      ,fHistPionEtaPt(0) 
						      ,fHistPionAccEtaPt(0) 
						      ,fHistPionRecEtaPt(0) 
						      ,fHistPionMatchEtaPt(0) 

						      ,fHistPiontotET(0) 
						      ,fHistPionAcctotET(0) 
						      ,fHistPionRectotET(0) 
						      ,fHistPionMatchtotET(0) 

						      ,fHistPionRectotETDep(0) 
						      ,fHistPionMatchtotETDep(0) 

						      ,fHistPionRecEtaEDepETDep(0) 
						      ,fHistPionMatchEtaEDepETDep(0) 

						      ,fHistPionRecEtaPtETDep(0) 
						      ,fHistPionMatchEtaPtETDep(0) 

						      ,fHistPionRecEtaETDep(0) 
						      ,fHistPionMatchEtaETDep(0) 

						      ,fHistPionRecResEET(0) 
						      ,fHistPionRecResPtET(0) 
						      ,fHistPionRecResE(0) 
						      ,fHistPionRecResPt(0) 
						      ,fHistPionRecResEDepETDep(0) 
						      ,fHistPionRecResPtETDep(0) 

						      ,fHistKaonEtaEET(0) 
						      ,fHistKaonAccEtaEET(0) 
						      ,fHistKaonRecEtaEET(0) 
						      ,fHistKaonMatchEtaEET(0) 

						      ,fHistKaonEtaPtET(0) 
						      ,fHistKaonAccEtaPtET(0) 
						      ,fHistKaonRecEtaPtET(0) 
						      ,fHistKaonMatchEtaPtET(0) 

						      ,fHistKaonEtaET(0) 
						      ,fHistKaonAccEtaET(0) 
						      ,fHistKaonRecEtaET(0) 
						      ,fHistKaonMatchEtaET(0) 

						      ,fHistKaonEtaE(0) 
						      ,fHistKaonAccEtaE(0) 
						      ,fHistKaonRecEtaE(0) 
						      ,fHistKaonMatchEtaE(0) 

						      ,fHistKaonEtaPt(0) 
						      ,fHistKaonAccEtaPt(0) 
						      ,fHistKaonRecEtaPt(0) 
						      ,fHistKaonMatchEtaPt(0) 

						      ,fHistKaontotET(0) 
						      ,fHistKaonAcctotET(0) 
						      ,fHistKaonRectotET(0) 
						      ,fHistKaonMatchtotET(0) 

						      ,fHistKaonRectotETDep(0) 
						      ,fHistKaonMatchtotETDep(0) 

						      ,fHistKaonRecEtaEDepETDep(0) 
						      ,fHistKaonMatchEtaEDepETDep(0) 

						      ,fHistKaonRecEtaPtETDep(0) 
						      ,fHistKaonMatchEtaPtETDep(0) 

						      ,fHistKaonRecEtaETDep(0) 
						      ,fHistKaonMatchEtaETDep(0) 

						      ,fHistKaonRecResEET(0) 
						      ,fHistKaonRecResPtET(0) 
						      ,fHistKaonRecResE(0) 
						      ,fHistKaonRecResPt(0) 

						      ,fHistKaonRecResEDepETDep(0) 
						      ,fHistKaonRecResPtETDep(0) 

						      ,fHistProtonEtaEET(0) 
						      ,fHistProtonAccEtaEET(0) 
						      ,fHistProtonRecEtaEET(0) 
						      ,fHistProtonMatchEtaEET(0) 

						      ,fHistProtonEtaPtET(0) 
						      ,fHistProtonAccEtaPtET(0) 
						      ,fHistProtonRecEtaPtET(0) 
						      ,fHistProtonMatchEtaPtET(0) 

						      ,fHistProtonEtaET(0) 
						      ,fHistProtonAccEtaET(0) 
						      ,fHistProtonRecEtaET(0) 
						      ,fHistProtonMatchEtaET(0) 

						      ,fHistProtonEtaE(0) 
						      ,fHistProtonAccEtaE(0) 
						      ,fHistProtonRecEtaE(0) 
						      ,fHistProtonMatchEtaE(0) 

						      ,fHistProtonEtaPt(0) 
						      ,fHistProtonAccEtaPt(0) 
						      ,fHistProtonRecEtaPt(0) 
						      ,fHistProtonMatchEtaPt(0) 

						      ,fHistProtontotET(0) 
						      ,fHistProtonAcctotET(0) 
						      ,fHistProtonRectotET(0) 
						      ,fHistProtonMatchtotET(0) 

						      ,fHistProtonRectotETDep(0) 
						      ,fHistProtonMatchtotETDep(0) 

						      ,fHistProtonRecEtaEDepETDep(0) 
						      ,fHistProtonMatchEtaEDepETDep(0) 

						      ,fHistProtonRecEtaPtETDep(0) 
						      ,fHistProtonMatchEtaPtETDep(0) 

						      ,fHistProtonRecEtaETDep(0) 
						      ,fHistProtonMatchEtaETDep(0) 

						      ,fHistProtonRecResEET(0) 
						      ,fHistProtonRecResPtET(0) 
						      ,fHistProtonRecResE(0) 
						      ,fHistProtonRecResPt(0) 

						      ,fHistProtonRecResEDepETDep(0) 
						      ,fHistProtonRecResPtETDep(0) 

						      ,fHistTotChargedtotET(0)
						      ,fHistTotChargedAcctotET(0)
						      ,fHistTotChargedRectotET(0)
						      ,fHistTotChargedRectotETDep(0)
						      ,fHistTotChargedMatchtotET(0)
						      ,fHistTotChargedMatchtotETDep(0)

						      ,fHistNeutronEtaEET(0) 
						      ,fHistNeutronAccEtaEET(0) 
						      ,fHistNeutronRecEtaEET(0) 

						      ,fHistNeutronEtaPtET(0) 
						      ,fHistNeutronAccEtaPtET(0) 
						      ,fHistNeutronRecEtaPtET(0) 

						      ,fHistNeutronEtaET(0) 
						      ,fHistNeutronAccEtaET(0) 
						      ,fHistNeutronRecEtaET(0) 

						      ,fHistNeutronEtaE(0) 
						      ,fHistNeutronAccEtaE(0) 
						      ,fHistNeutronRecEtaE(0) 

						      ,fHistNeutronEtaPt(0) 
						      ,fHistNeutronAccEtaPt(0) 
						      ,fHistNeutronRecEtaPt(0) 

						      ,fHistNeutrontotET(0) 
						      ,fHistNeutronAcctotET(0) 
						      ,fHistNeutronRectotET(0) 

						      ,fHistNeutronRectotETDep(0)

						      ,fHistNeutronRecEtaEDepETDep(0) 
						      ,fHistNeutronRecEtaETDep(0) 
						      ,fHistNeutronRecEtaPtETDep(0) 

						      ,fHistK0EtaEET(0) 
						      ,fHistK0RecEtaEET(0) 

						      ,fHistK0EtaPtET(0) 
						      ,fHistK0RecEtaPtET(0) 

						      ,fHistK0EtaET(0) 
						      ,fHistK0RecEtaET(0) 

						      ,fHistK0EtaE(0) 
						      ,fHistK0RecEtaE(0) 

						      ,fHistK0EtaPt(0) 
						      ,fHistK0RecEtaPt(0) 

						      ,fHistK0totET(0) 
						      ,fHistK0RectotET(0) 
						      ,fHistK0RectotETDep(0) 

						      ,fHistK0RecEtaEDepETDep(0) 
						      ,fHistK0RecEtaETDep(0) 

						      ,fHistK0RecEtaPtETDep(0) 

						      ,fHistLambdaEtaEET(0) 
						      ,fHistLambdaRecEtaEET(0) 

						      ,fHistLambdaEtaPtET(0) 
						      ,fHistLambdaRecEtaPtET(0) 

						      ,fHistLambdaEtaET(0) 
						      ,fHistLambdaRecEtaET(0) 

						      ,fHistLambdaEtaE(0) 
						      ,fHistLambdaRecEtaE(0) 

						      ,fHistLambdaEtaPt(0) 
						      ,fHistLambdaRecEtaPt(0) 

						      ,fHistLambdatotET(0) 
						      ,fHistLambdaRectotET(0) 
						      ,fHistLambdaRectotETDep(0) 

						      ,fHistLambdaRecEtaEDepETDep(0) 
						      ,fHistLambdaRecEtaETDep(0) 

						      ,fHistLambdaRecEtaPtETDep(0) 

						      ,fHistTotNeutraltotET(0)
						      ,fHistTotNeutralRectotET(0)
						      ,fHistTotNeutralRectotETDep(0)

						      ,fHistTotaltotET(0)
						      ,fHistTotalAcctotET(0)
						      ,fHistTotalRectotET(0)
						      ,fHistTotalRectotETDep(0)

						      ,fHistElectronFirstMother(0) 
						      ,fHistElectronFirstMotherXY(0) 
						      ,fHistElectronNDaughters(0) 
						      ,fHistElectronDaughters(0) 
						      ,fHistElectronDaughtersXY(0) 

						      ,fHistElectronFirstMotherAcc(0)  
						      ,fHistElectronFirstMotherXYAcc(0)  
						      ,fHistElectronNDaughtersAcc(0) 
						      ,fHistElectronDaughtersAcc(0) 
						      ,fHistElectronDaughtersXYAcc(0) 

						      ,fHistElectronFirstMotherRec(0)  
						      ,fHistElectronFirstMotherXYRec(0)  
						      ,fHistElectronNDaughtersRec(0) 
						      ,fHistElectronDaughtersRec(0) 
						      ,fHistElectronDaughtersXYRec(0) 

						      ,fHistNPPElectronFirstMother(0) 
						      ,fHistNPPElectronFirstMotherXY(0) 
						      ,fHistNPPElectronNDaughters(0) 
						      ,fHistNPPElectronDaughters(0) 
						      ,fHistNPPElectronDaughtersXY(0) 

						      ,fHistNPPElectronFirstMotherAcc(0)  
						      ,fHistNPPElectronFirstMotherXYAcc(0)  
						      ,fHistNPPElectronNDaughtersAcc(0) 
						      ,fHistNPPElectronDaughtersAcc(0) 
						      ,fHistNPPElectronDaughtersXYAcc(0) 

						      ,fHistNPPElectronFirstMotherRec(0)  
						      ,fHistNPPElectronFirstMotherXYRec(0)  
						      ,fHistNPPElectronNDaughtersRec(0) 
						      ,fHistNPPElectronDaughtersRec(0) 
						      ,fHistNPPElectronDaughtersXYRec(0) 

						      ,fHistGammaFirstMother(0) 
						      ,fHistGammaFirstMotherXY(0) 
						      ,fHistGammaNDaughters(0) 
						      ,fHistGammaDaughters(0) 
						      ,fHistGammaDaughtersXY(0) 
						      ,fHistConvGammaDaughtersXY(0) 
						      ,fHistNonConvGammaDaughtersXY(0) 

						      ,fHistGammaFirstMotherAcc(0)  
						      ,fHistGammaFirstMotherXYAcc(0)  
						      ,fHistGammaNDaughtersAcc(0) 
						      ,fHistGammaDaughtersAcc(0) 
						      ,fHistGammaDaughtersXYAcc(0) 
						      ,fHistConvGammaDaughtersXYAcc(0) 
						      ,fHistNonConvGammaDaughtersXYAcc(0) 

						      ,fHistGammaFirstMotherRec(0)  
						      ,fHistGammaFirstMotherXYRec(0)  
						      ,fHistGammaNDaughtersRec(0) 
						      ,fHistGammaDaughtersRec(0) 
						      ,fHistGammaDaughtersXYRec(0) 
						      ,fHistConvGammaDaughtersXYRec(0) 
						      ,fHistNonConvGammaDaughtersXYRec(0) 

						      ,fHistNPPGammaFirstMother(0) 
						      ,fHistNPPGammaFirstMotherXY(0) 
						      ,fHistNPPGammaNDaughters(0) 
						      ,fHistNPPGammaDaughters(0) 
						      ,fHistNPPGammaDaughtersXY(0) 

						      ,fHistNPPGammaFirstMotherAcc(0)  
						      ,fHistNPPGammaFirstMotherXYAcc(0)  
						      ,fHistNPPGammaNDaughtersAcc(0) 
						      ,fHistNPPGammaDaughtersAcc(0) 
						      ,fHistNPPGammaDaughtersXYAcc(0) 

						      ,fHistNPPGammaFirstMotherRec(0)  
						      ,fHistNPPGammaFirstMotherXYRec(0)  
						      ,fHistNPPGammaNDaughtersRec(0) 
						      ,fHistNPPGammaDaughtersRec(0) 
						      ,fHistNPPGammaDaughtersXYRec(0) 

						      ,fHistAllERecEMC(0)	
						      ,fHistAllPtRecPtMC(0)
						      ,fHistElectronERecEMC(0)	
						      ,fHistGammaERecEMC(0)

						      ,fHistChargedRes(0)
						      ,fHistChargedRes2(0)
						      ,fHistChargedRes3(0)
						      ,fHistNeutralRes(0)
						      ,fHistElectronRes(0)
						      ,fHistGammaRes(0)

						      ,fHistIsInAcc(0)
{//constructor
  fHistogramNameSuffix = TString("EmcalMC");
	
  fResCut = 0.02;
  //fResCut = fEmcalTrackDistanceCut;
	
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));
  //TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));
  TGeoManager::Import("geometry.root");
  //fGeoUt = new AliEMCALGeometry("EMCAL_FIRSTYEAR","EMCAL");
}

// dtor
AliAnalysisEmEtMonteCarlo::~AliAnalysisEmEtMonteCarlo() 
{//destructor
  //Marcelo, have you really cleaned up all memory here?  What about the histos?
  delete fGeoUt;
  delete fHistPrimEtaEET;
  delete fHistPrimEtaPtET; 
  delete fHistPrimEtaET;
  delete fHistPrimtotET;
	
  delete fHistPrimAccEtaEET; 
  delete fHistPrimAccEtaPtET; 
  delete fHistPrimAccEtaET; 
  delete fHistPrimAcctotET;
	
  delete fHistPrimRecEtaEET; 
  delete fHistPrimRecEtaPtET; 
  delete fHistPrimRecEtaET; 
  delete fHistPrimRectotET;

  delete fHistPrimRecEtaEDepETDep; 
  delete fHistPrimRecEtaPtETDep; 
  delete fHistPrimRecEtaETDep; 
  delete fHistPrimRectotETDep;
	
  // *******************
  // electron ET
  // *******************
  delete fHistElectronEtaEET; 
  delete fHistElectronEtaPtET; 
  delete fHistElectronEtaET; 
  delete fHistElectronEtaE; 
  delete fHistElectronEtaPt; 
  delete fHistElectrontotET; 

  delete fHistConvElectronEtaEET;  
  delete fHistConvElectronEtaPtET;  
  delete fHistConvElectronEtaET;  
  delete fHistConvElectronEtaE;  
  delete fHistConvElectronEtaPt;  
  delete fHistConvElectrontotET;  

  delete fHistScatElectronEtaEET;  
  delete fHistScatElectronEtaPtET;  
  delete fHistScatElectronEtaET;  
  delete fHistScatElectronEtaE;  
  delete fHistScatElectronEtaPt;  
  delete fHistScatElectrontotET;  
	
  // *******************
  // total electron ET
  // *******************
  delete fHistTotElectrontotET;
	
  // *******************
  // gamma ET
  // *******************
  delete fHistGammaEtaEET;  
  delete fHistGammaEtaPtET;  
  delete fHistGammaEtaET;  
  delete fHistGammaEtaE;  
  delete fHistGammaEtaPt;  
  delete fHistGammatotET;  
	
  delete fHistAnnihGammaEtaEET;  
  delete fHistAnnihGammaEtaPtET;  
  delete fHistAnnihGammaEtaET;  
  delete fHistAnnihGammaEtaE;  
  delete fHistAnnihGammaEtaPt;  
  delete fHistAnnihGammatotET;  

  delete fHistScatGammaEtaEET;  
  delete fHistScatGammaEtaPtET;  
  delete fHistScatGammaEtaET;  
  delete fHistScatGammaEtaE;  
  delete fHistScatGammaEtaPt;  
  delete fHistScatGammatotET;  

  delete fHistConvGammaEtaEET;  
  delete fHistConvGammaEtaPtET;  
  delete fHistConvGammaEtaET;  
  delete fHistConvGammaEtaE;  
  delete fHistConvGammaEtaPt;  
  delete fHistConvGammatotET;  
	
  delete fHistNonConvGammaEtaEET;  
  delete fHistNonConvGammaEtaPtET;  
  delete fHistNonConvGammaEtaET;  
  delete fHistNonConvGammaEtaE;  
  delete fHistNonConvGammaEtaPt;  
  delete fHistNonConvGammatotET;  
	
  // *******************
  // total gamma ET
  // *******************
  delete fHistTotGammatotET;

  // *******************
  // total electromagnetic ET
  // *******************
  delete fHistTotEMtotET;

  // non-primary electromagnetic ET
  delete fHistNPPElectronEtaEET; 
  delete fHistNPPElectronEtaPtET; 
  delete fHistNPPElectronEtaET; 
  delete fHistNPPElectronEtaE; 
  delete fHistNPPElectronEtaPt; 
  delete fHistNPPElectrontotET; 

  delete fHistNPPGammaEtaEET; 
  delete fHistNPPGammaEtaPtET; 
  delete fHistNPPGammaEtaET; 
  delete fHistNPPGammaEtaE; 
  delete fHistNPPGammaEtaPt; 
  delete fHistNPPGammatotET; 

  delete fHistTotNPPEMtotET;

  delete fHistNPPPi0GammaEtaEET; 
  delete fHistNPPPi0GammaEtaPtET; 
  delete fHistNPPPi0GammaEtaET; 
  delete fHistNPPPi0GammaEtaE; 
  delete fHistNPPPi0GammaEtaPt; 
  delete fHistNPPPi0GammatotET; 
		
  // *******************
  // electron ET inside EMCal acceptance
  // *******************
  delete fHistElectronAccEtaEET; 
  delete fHistElectronAccEtaPtET; 
  delete fHistElectronAccEtaET; 
  delete fHistElectronAccEtaE; 
  delete fHistElectronAccEtaPt; 
  delete fHistElectronAcctotET; 
	
  delete fHistConvElectronAccEtaEET;  
  delete fHistConvElectronAccEtaPtET;  
  delete fHistConvElectronAccEtaET;  
  delete fHistConvElectronAccEtaE;  
  delete fHistConvElectronAccEtaPt;  
  delete fHistConvElectronAcctotET;  
	
  delete fHistScatElectronAccEtaEET;  
  delete fHistScatElectronAccEtaPtET;  
  delete fHistScatElectronAccEtaET;  
  delete fHistScatElectronAccEtaE;  
  delete fHistScatElectronAccEtaPt;  
  delete fHistScatElectronAcctotET;  
	
  // *******************
  // total electron ET inside EMCal acceptance
  // *******************
  delete fHistTotElectronAcctotET;

  // *******************
  // gamma ET inside EMCal acceptance
  // *******************
  delete fHistGammaAccEtaEET;  
  delete fHistGammaAccEtaPtET;  
  delete fHistGammaAccEtaET;  
  delete fHistGammaAccEtaE;  
  delete fHistGammaAccEtaPt;  
  delete fHistGammaAcctotET;  
	
  delete fHistAnnihGammaAccEtaEET;  
  delete fHistAnnihGammaAccEtaPtET;  
  delete fHistAnnihGammaAccEtaET;  
  delete fHistAnnihGammaAccEtaE;  
  delete fHistAnnihGammaAccEtaPt;  
  delete fHistAnnihGammaAcctotET;  
	
  delete fHistScatGammaAccEtaEET;  
  delete fHistScatGammaAccEtaPtET;  
  delete fHistScatGammaAccEtaET;  
  delete fHistScatGammaAccEtaE;  
  delete fHistScatGammaAccEtaPt;  
  delete fHistScatGammaAcctotET;  
	
  delete fHistConvGammaAccEtaEET;  
  delete fHistConvGammaAccEtaPtET;  
  delete fHistConvGammaAccEtaET;  
  delete fHistConvGammaAccEtaE;  
  delete fHistConvGammaAccEtaPt;  
  delete fHistConvGammaAcctotET;  
	
  delete fHistNonConvGammaAccEtaEET;  
  delete fHistNonConvGammaAccEtaPtET;  
  delete fHistNonConvGammaAccEtaET;  
  delete fHistNonConvGammaAccEtaE;  
  delete fHistNonConvGammaAccEtaPt;  
  delete fHistNonConvGammaAcctotET;  
	
  // *******************
  // total gamma ET inside EMCal acceptance
  // *******************
  delete fHistTotGammaAcctotET;

  // *******************
  // total electromagnetic ET inside EMCal acceptance
  // *******************
  delete fHistTotEMAcctotET;

  // non-primary electromagnetic ET
  delete fHistNPPElectronAccEtaEET; 
  delete fHistNPPElectronAccEtaPtET; 
  delete fHistNPPElectronAccEtaE; 
  delete fHistNPPElectronAccEtaPt; 
	
  delete fHistNPPGammaAccEtaEET; 
  delete fHistNPPGammaAccEtaPtET; 
  delete fHistNPPGammaAccEtaE; 
  delete fHistNPPGammaAccEtaPt; 	
	
  // *******************
  // electron ET reconstructed in EMCal
  // *******************
  delete fHistElectronRecEtaEET; 
  delete fHistElectronRecEtaPtET; 
  delete fHistElectronRecEtaET; 
  delete fHistElectronRecEtaE; 
  delete fHistElectronRecEtaPt; 
  delete fHistElectronRectotET; 
	
  delete fHistConvElectronRecEtaEET;  
  delete fHistConvElectronRecEtaPtET;  
  delete fHistConvElectronRecEtaET;  
  delete fHistConvElectronRecEtaE;  
  delete fHistConvElectronRecEtaPt;  
  delete fHistConvElectronRectotET;  
	
  delete fHistScatElectronRecEtaEET;  
  delete fHistScatElectronRecEtaPtET;  
  delete fHistScatElectronRecEtaET;  
  delete fHistScatElectronRecEtaE;  
  delete fHistScatElectronRecEtaPt;  
  delete fHistScatElectronRectotET;  
	
  // *******************
  // total Electron ET reconstructed in EMCal
  // *******************
  delete fHistTotElectronRectotET;

  // *******************
  // gamma ET reconstructed in EMCal
  // *******************
  delete fHistGammaRecEtaEET;  
  delete fHistGammaRecEtaPtET;  
  delete fHistGammaRecEtaET;  
  delete fHistGammaRecEtaE;  
  delete fHistGammaRecEtaPt;  
  delete fHistGammaRectotET;  
	
  delete fHistAnnihGammaRecEtaEET;  
  delete fHistAnnihGammaRecEtaPtET;  
  delete fHistAnnihGammaRecEtaET;  
  delete fHistAnnihGammaRecEtaE;  
  delete fHistAnnihGammaRecEtaPt;  
  delete fHistAnnihGammaRectotET;  
	
  delete fHistScatGammaRecEtaEET;  
  delete fHistScatGammaRecEtaPtET;  
  delete fHistScatGammaRecEtaET;  
  delete fHistScatGammaRecEtaE;  
  delete fHistScatGammaRecEtaPt;  
  delete fHistScatGammaRectotET;  

  // *******************
  // total gamma ET reconstructed in EMCal
  // *******************
  delete fHistTotGammaRectotET;

  // *******************
  // total EM ET reconstructed in EMCal
  // *******************
  delete fHistTotEMRectotET;

  // non-primary electromagnetic ET
  delete fHistNPPElectronRecEtaEET; 
  delete fHistNPPElectronRecEtaPtET; 
  delete fHistNPPElectronRecEtaET; 
  delete fHistNPPElectronRecEtaE; 
  delete fHistNPPElectronRecEtaPt; 
  delete fHistNPPElectronRectotET; 
	
  delete fHistNPPGammaRecEtaEET; 
  delete fHistNPPGammaRecEtaPtET; 
  delete fHistNPPGammaRecEtaET; 
  delete fHistNPPGammaRecEtaE; 
  delete fHistNPPGammaRecEtaPt; 
  delete fHistNPPGammaRectotET; 
	
  delete fHistTotNPPEMRectotET;

  delete fHistNPPPi0GammaRecEtaEET; 
  delete fHistNPPPi0GammaRecEtaPtET; 
  delete fHistNPPPi0GammaRecEtaET; 
  delete fHistNPPPi0GammaRecEtaE; 
  delete fHistNPPPi0GammaRecEtaPt; 
  delete fHistNPPPi0GammaRectotET; 
	
  // *******************
  // muon ET (+ and -)
  // *******************
  delete fHistMuonEtaEET; 
  delete fHistMuonAccEtaEET; 
  delete fHistMuonRecEtaEET; 
  delete fHistMuonMatchEtaEET; 

  delete fHistMuonEtaPtET; 
  delete fHistMuonAccEtaPtET; 
  delete fHistMuonRecEtaPtET; 
  delete fHistMuonMatchEtaPtET; 

  delete fHistMuonEtaET; 
  delete fHistMuonAccEtaET; 
  delete fHistMuonRecEtaET; 
  delete fHistMuonMatchEtaET; 
	
  delete fHistMuonEtaE; 
  delete fHistMuonAccEtaE; 
  delete fHistMuonRecEtaE; 
  delete fHistMuonMatchEtaE; 
	
  delete fHistMuonEtaPt; 
  delete fHistMuonAccEtaPt; 
  delete fHistMuonRecEtaPt; 
  delete fHistMuonMatchEtaPt; 
	
  delete fHistMuontotET; 
  delete fHistMuonAcctotET; 
  delete fHistMuonRectotET; 
  delete fHistMuonMatchtotET; 
	
  delete fHistMuonRectotETDep; 
  delete fHistMuonMatchtotETDep; 
	
  delete fHistMuonRecEtaEDepETDep; 
  delete fHistMuonMatchEtaEDepETDep; 

  delete fHistMuonRecEtaPtETDep; 
  delete fHistMuonMatchEtaPtETDep; 
	
  delete fHistMuonRecEtaETDep; 
  delete fHistMuonMatchEtaETDep; 

  delete fHistMuonRecResEET; 
  delete fHistMuonRecResPtET; 
  delete fHistMuonRecResE; 
  delete fHistMuonRecResPt; 
  delete fHistMuonRecResEDepETDep; 
  delete fHistMuonRecResPtETDep; 
	
  // *******************
  // pion ET (+ and -)
  // *******************
  delete fHistPionEtaEET; 
  delete fHistPionAccEtaEET; 
  delete fHistPionRecEtaEET; 
  delete fHistPionMatchEtaEET; 
	
  delete fHistPionEtaPtET; 
  delete fHistPionAccEtaPtET; 
  delete fHistPionRecEtaPtET; 
  delete fHistPionMatchEtaPtET; 
	
  delete fHistPionEtaET; 
  delete fHistPionAccEtaET; 
  delete fHistPionRecEtaET; 
  delete fHistPionMatchEtaET; 
	
  delete fHistPionEtaE; 
  delete fHistPionAccEtaE; 
  delete fHistPionRecEtaE; 
  delete fHistPionMatchEtaE; 
	
  delete fHistPionEtaPt; 
  delete fHistPionAccEtaPt; 
  delete fHistPionRecEtaPt; 
  delete fHistPionMatchEtaPt; 
	
  delete fHistPiontotET; 
  delete fHistPionAcctotET; 
  delete fHistPionRectotET; 
  delete fHistPionMatchtotET; 
	
  delete fHistPionRectotETDep; 
  delete fHistPionMatchtotETDep; 
	
  delete fHistPionRecEtaEDepETDep; 
  delete fHistPionMatchEtaEDepETDep; 

  delete fHistPionRecEtaPtETDep; 
  delete fHistPionMatchEtaPtETDep; 
	
  delete fHistPionRecEtaETDep; 
  delete fHistPionMatchEtaETDep; 
	
  delete fHistPionRecResEET; 
  delete fHistPionRecResPtET; 
  delete fHistPionRecResE; 
  delete fHistPionRecResPt; 
  delete fHistPionRecResEDepETDep; 
  delete fHistPionRecResPtETDep; 
	
  // *******************
  // charged kaon (+ and -) ET
  // *******************
  delete fHistKaonEtaEET; 
  delete fHistKaonAccEtaEET; 
  delete fHistKaonRecEtaEET; 
  delete fHistKaonMatchEtaEET; 
	
  delete fHistKaonEtaPtET; 
  delete fHistKaonAccEtaPtET; 
  delete fHistKaonRecEtaPtET; 
  delete fHistKaonMatchEtaPtET; 
	
  delete fHistKaonEtaET; 
  delete fHistKaonAccEtaET; 
  delete fHistKaonRecEtaET; 
  delete fHistKaonMatchEtaET; 
	
  delete fHistKaonEtaE; 
  delete fHistKaonAccEtaE; 
  delete fHistKaonRecEtaE; 
  delete fHistKaonMatchEtaE; 
	
  delete fHistKaonEtaPt; 
  delete fHistKaonAccEtaPt; 
  delete fHistKaonRecEtaPt; 
  delete fHistKaonMatchEtaPt; 

  delete fHistKaontotET; 
  delete fHistKaonAcctotET; 
  delete fHistKaonRectotET; 
  delete fHistKaonMatchtotET; 
	
  delete fHistKaonRectotETDep; 
  delete fHistKaonMatchtotETDep; 
	
  delete fHistKaonRecEtaEDepETDep; 
  delete fHistKaonMatchEtaEDepETDep; 

  delete fHistKaonRecEtaPtETDep; 
  delete fHistKaonMatchEtaPtETDep; 
	
  delete fHistKaonRecEtaETDep; 
  delete fHistKaonMatchEtaETDep; 
	
  delete fHistKaonRecResEET; 
  delete fHistKaonRecResPtET; 
  delete fHistKaonRecResE; 
  delete fHistKaonRecResPt; 
  delete fHistKaonRecResEDepETDep; 
  delete fHistKaonRecResPtETDep; 	
	
  // *******************
  // proton (anti) ET
  // *******************
  delete fHistProtonEtaEET; 
  delete fHistProtonAccEtaEET; 
  delete fHistProtonRecEtaEET; 
  delete fHistProtonMatchEtaEET; 
	
  delete fHistProtonEtaPtET; 
  delete fHistProtonAccEtaPtET; 
  delete fHistProtonRecEtaPtET; 
  delete fHistProtonMatchEtaPtET; 
	
  delete fHistProtonEtaET; 
  delete fHistProtonAccEtaET; 
  delete fHistProtonRecEtaET; 
  delete fHistProtonMatchEtaET; 
	
  delete fHistProtonEtaE; 
  delete fHistProtonAccEtaE; 
  delete fHistProtonRecEtaE; 
  delete fHistProtonMatchEtaE; 
	
  delete fHistProtonEtaPt; 
  delete fHistProtonAccEtaPt; 
  delete fHistProtonRecEtaPt; 
  delete fHistProtonMatchEtaPt; 

  delete fHistProtontotET; 
  delete fHistProtonAcctotET; 
  delete fHistProtonRectotET; 
  delete fHistProtonMatchtotET; 
	
  delete fHistProtonRectotETDep; 
  delete fHistProtonMatchtotETDep; 
	
  delete fHistProtonRecEtaEDepETDep; 
  delete fHistProtonMatchEtaEDepETDep; 
	
  delete fHistProtonRecEtaPtETDep; 
  delete fHistProtonMatchEtaPtETDep; 
	
  delete fHistProtonRecEtaETDep; 
  delete fHistProtonMatchEtaETDep; 

  delete fHistProtonRecResEET; 
  delete fHistProtonRecResPtET; 
  delete fHistProtonRecResE; 
  delete fHistProtonRecResPt; 
  delete fHistProtonRecResEDepETDep; 
  delete fHistProtonRecResPtETDep; 
	
  // *******************
  // total charged ET
  // *******************
  delete fHistTotChargedtotET;
  delete fHistTotChargedAcctotET;
  delete fHistTotChargedRectotET;
  delete fHistTotChargedRectotETDep;
  delete fHistTotChargedMatchtotET;
  delete fHistTotChargedMatchtotETDep;
	
  // *******************
  // neutron (anti) ET
  // *******************
  delete fHistNeutronEtaEET; 
  delete fHistNeutronAccEtaEET; 
  delete fHistNeutronRecEtaEET; 
	
  delete fHistNeutronEtaPtET; 
  delete fHistNeutronAccEtaPtET; 
  delete fHistNeutronRecEtaPtET; 
	
  delete fHistNeutronEtaET; 
  delete fHistNeutronAccEtaET; 
  delete fHistNeutronRecEtaET; 
	
  delete fHistNeutronEtaE; 
  delete fHistNeutronAccEtaE; 
  delete fHistNeutronRecEtaE; 
	
  delete fHistNeutronEtaPt; 
  delete fHistNeutronAccEtaPt; 
  delete fHistNeutronRecEtaPt; 
	
  delete fHistNeutrontotET; 
  delete fHistNeutronAcctotET; 
  delete fHistNeutronRectotET; 
  delete fHistNeutronRectotETDep; 
	
  delete fHistNeutronRecEtaEDepETDep; 
  delete fHistNeutronRecEtaETDep; 
	
  delete fHistNeutronRecEtaPtETDep; 
		
  // *******************
  // neutral kaon ET
  // *******************
  delete fHistK0EtaEET; 
  delete fHistK0RecEtaEET; 
	
  delete fHistK0EtaPtET; 
  delete fHistK0RecEtaPtET; 
	
  delete fHistK0EtaET; 
  delete fHistK0RecEtaET; 
	
  delete fHistK0EtaE; 
  delete fHistK0RecEtaE; 
	
  delete fHistK0EtaPt; 
  delete fHistK0RecEtaPt; 

  delete fHistK0totET; 
  delete fHistK0RectotET; 
	
  delete fHistK0RectotETDep; 
	
  delete fHistK0RecEtaEDepETDep; 
  delete fHistK0RecEtaETDep; 
	
  delete fHistK0RecEtaPtETDep; 
		
  // *******************
  // Lambda(anti) ET
  // *******************
  delete fHistLambdaEtaEET; 
  delete fHistLambdaRecEtaEET; 
	
  delete fHistLambdaEtaPtET; 
  delete fHistLambdaRecEtaPtET; 
	
  delete fHistLambdaEtaET; 
  delete fHistLambdaRecEtaET; 
	
  delete fHistLambdaEtaE; 
  delete fHistLambdaRecEtaE; 
	
  delete fHistLambdaEtaPt; 
  delete fHistLambdaRecEtaPt; 
	
  delete fHistLambdatotET; 
  delete fHistLambdaRectotET; 
	
  delete fHistLambdaRectotETDep; 
	
  delete fHistLambdaRecEtaEDepETDep; 
  delete fHistLambdaRecEtaETDep; 
	
  delete fHistLambdaRecEtaPtETDep; 

  // *******************
  // total neutral ET
  // *******************
  delete fHistTotNeutraltotET;
  delete fHistTotNeutralRectotET;
  delete fHistTotNeutralRectotETDep;
	
  // *******************
  // total ET
  // *******************
  delete fHistTotaltotET;
  delete fHistTotalAcctotET;
  delete fHistTotalRectotET;
  delete fHistTotalRectotETDep;
	
  // *******************
  // some checks
  // *******************

  // check produced electrons
  delete fHistElectronFirstMother; 
  delete fHistElectronFirstMotherXY; 
  delete fHistElectronNDaughters; 
  delete fHistElectronDaughters; 
  delete fHistElectronDaughtersXY; 

  delete fHistElectronFirstMotherAcc;  
  delete fHistElectronFirstMotherXYAcc;  
  delete fHistElectronNDaughtersAcc; 
  delete fHistElectronDaughtersAcc; 
  delete fHistElectronDaughtersXYAcc; 

  delete fHistElectronFirstMotherRec;  
  delete fHistElectronFirstMotherXYRec;  
  delete fHistElectronNDaughtersRec; 
  delete fHistElectronDaughtersRec; 
  delete fHistElectronDaughtersXYRec; 

  delete fHistNPPElectronFirstMother; 
  delete fHistNPPElectronFirstMotherXY; 
  delete fHistNPPElectronNDaughters; 
  delete fHistNPPElectronDaughters; 
  delete fHistNPPElectronDaughtersXY; 
	
  delete fHistNPPElectronFirstMotherAcc;  
  delete fHistNPPElectronFirstMotherXYAcc;  
  delete fHistNPPElectronNDaughtersAcc; 
  delete fHistNPPElectronDaughtersAcc; 
  delete fHistNPPElectronDaughtersXYAcc; 
	
  delete fHistNPPElectronFirstMotherRec;  
  delete fHistNPPElectronFirstMotherXYRec;  
  delete fHistNPPElectronNDaughtersRec; 
  delete fHistNPPElectronDaughtersRec; 
  delete fHistNPPElectronDaughtersXYRec; 
	
  // check produced gammas
  delete fHistGammaFirstMother; 
  delete fHistGammaFirstMotherXY; 
  delete fHistGammaNDaughters; 
  delete fHistGammaDaughters; 
  delete fHistGammaDaughtersXY; 
  delete fHistConvGammaDaughtersXY; 
  delete fHistNonConvGammaDaughtersXY; 
	
  delete fHistGammaFirstMotherAcc;  
  delete fHistGammaFirstMotherXYAcc;  
  delete fHistGammaNDaughtersAcc; 
  delete fHistGammaDaughtersAcc; 
  delete fHistGammaDaughtersXYAcc; 
  delete fHistConvGammaDaughtersXYAcc; 
  delete fHistNonConvGammaDaughtersXYAcc; 
	
  delete fHistGammaFirstMotherRec;  
  delete fHistGammaFirstMotherXYRec;  
  delete fHistGammaNDaughtersRec; 
  delete fHistGammaDaughtersRec; 
  delete fHistGammaDaughtersXYRec; 
  delete fHistConvGammaDaughtersXYRec; 
  delete fHistNonConvGammaDaughtersXYRec; 
	
  delete fHistNPPGammaFirstMother; 
  delete fHistNPPGammaFirstMotherXY; 
  delete fHistNPPGammaNDaughters; 
  delete fHistNPPGammaDaughters; 
  delete fHistNPPGammaDaughtersXY; 
	
  delete fHistNPPGammaFirstMotherAcc;  
  delete fHistNPPGammaFirstMotherXYAcc;  
  delete fHistNPPGammaNDaughtersAcc; 
  delete fHistNPPGammaDaughtersAcc; 
  delete fHistNPPGammaDaughtersXYAcc; 
	
  delete fHistNPPGammaFirstMotherRec;  
  delete fHistNPPGammaFirstMotherXYRec;  
  delete fHistNPPGammaNDaughtersRec; 
  delete fHistNPPGammaDaughtersRec; 
  delete fHistNPPGammaDaughtersXYRec; 

  //check projections
  delete fHistAllERecEMC;	
  delete fHistAllPtRecPtMC;
  delete fHistElectronERecEMC;	
  delete fHistGammaERecEMC;
	
  delete fHistChargedRes;
  delete fHistChargedRes2;
  delete fHistChargedRes3;
  delete fHistNeutralRes;
  delete fHistElectronRes;
  delete fHistGammaRes;
	
  delete fHistIsInAcc;
  //delete   TH2F * yyyyyy
}

Int_t AliAnalysisEmEtMonteCarlo::AnalyseEvent(AliVEvent* ev)
{ // analyse MC event
  //ResetEventValues();
	
  // Get us an mc event
  if(!ev)
    {
      Printf("ERROR: Event does not exist");   
      return 0;
    }
  AliMCEvent *event = dynamic_cast<AliMCEvent*>(ev);
		
  // Hijing header
  AliGenEventHeader* genHeader = event->GenEventHeader();
  if(!genHeader){
    Printf("ERROR: Event generation header does not exist");   
    return 0;
  }
    	
  // Let's play with the stack!
  AliStack *stack = event->Stack();

  if (!stack)
    {
      Printf("ERROR: Could not get stack");
      return 0;
    }	
	
  //Int_t nStackTracks = stack->GetNtrack();
  Int_t nStackTracks = event->GetNumberOfTracks();
	
  for (Int_t iPart = 0; iPart < nStackTracks; iPart++){
    AliMCParticle* aliPart = (AliMCParticle*)event->GetTrack(iPart);

    //TParticle *part = stack->Particle(iPart);
    TParticle *part = aliPart->Particle();
    TParticle *partMom = 0;
    TParticle *partDaughter = 0;
		
    if (!part){
      Printf("ERROR: Could not get particle %d", iPart);
      continue;
    }

    Int_t iPartMom = part->GetMother(0);
    Int_t iPartDaughter = 0;
    Int_t nPartDaughters = part->GetNDaughters();
	
    TParticlePDG *pdg = part->GetPDG(0);
    TParticlePDG *pdgMom = 0;
    TParticlePDG *pdgDaugther = 0;
	
    if (!pdg){
      Printf("ERROR-1: Could not get particle PDG %d", iPart);
      continue;
    }		
		
    //create an external track param for projection
    AliExternalTrackParam* extParam = CreateExternalTrackParam(part);
	
    if ((iPartMom>=0) && (iPartMom < nStackTracks))
      {
	partMom = stack->Particle(iPartMom);
	pdgMom = partMom->GetPDG(0);
      }
		
    // Check if it is a primary particle
		
    // Check for reasonable (for now neutral and singly charged) charge on the particle
    //TODO:Maybe not only singly charged?
    if (TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloSingleChargedParticle())<1e-3 && TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloNeutralParticle())<1e-3) continue;
		
    fMultiplicity++;
		
    // Inside ALICE central barrel acceptance
		
    if (TMath::Abs(part->Eta()) < fCuts->GetCommonEtaCut())
      {
	Double_t et = CalcET(part,pdg);
			
	if (et < 0) continue;

	if (IsPrimary(stack,iPart,pdg,iPartMom,pdgMom))
	  {		
	    if (stack->IsPhysicalPrimary(iPart))
	      {
		fHistPrimEtaEET->Fill(part->Energy(),part->Eta(),et);
		fHistPrimEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		fHistPrimEtaET->Fill(et,part->Eta());											
		fPrimtotET += et;
					
		//if (IsInAcceptance(part,pdg,extParam)) 
		if (IsInAcceptance(aliPart)) 
		  {
		    fHistPrimAccEtaEET->Fill(part->Energy(),part->Eta(),et);
		    fHistPrimAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		    fHistPrimAccEtaET->Fill(et,part->Eta());																
		    fPrimAcctotET += et;
		  }				
	      }

	    // Fill up total E_T counters for each particle species
	    if (pdg->PdgCode() == fgProtonCode || pdg->PdgCode() == fgAntiProtonCode)
	      {
		fProtontotET += et;
		fHistProtonEtaEET->Fill(part->Energy(),part->Eta(),et);
		fHistProtonEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		fHistProtonEtaET->Fill(et,part->Eta());							
		fHistProtonEtaE->Fill(part->Energy(),part->Eta());
		fHistProtonEtaPt->Fill(part->Pt(),part->Eta());
					
		// inside EMCal acceptance
		//if (IsInAcceptance(part,pdg,extParam)) 
		if (IsInAcceptance(aliPart)) 
		  {
		    fProtonAcctotET += et;
		    fHistProtonAccEtaEET->Fill(part->Energy(),part->Eta(),et);
		    fHistProtonAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		    fHistProtonAccEtaET->Fill(et,part->Eta());
		    fHistProtonAccEtaE->Fill(part->Energy(),part->Eta());
		    fHistProtonAccEtaPt->Fill(part->Pt(),part->Eta());													
		  }						
	      }
	    if (pdg->PdgCode() == fgPiPlusCode || pdg->PdgCode() == fgPiMinusCode)
	      {
		fPiontotET += et;
		fHistPionEtaEET->Fill(part->Energy(),part->Eta(),et);
		fHistPionEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		fHistPionEtaET->Fill(et,part->Eta());
		fHistPionEtaE->Fill(part->Energy(),part->Eta());
		fHistPionEtaPt->Fill(part->Pt(),part->Eta());							
		// inside EMCal acceptance
		//if (IsInAcceptance(part,pdg,extParam)) 
		if (IsInAcceptance(aliPart)) 
		  {
		    fPionAcctotET += et;
		    fHistPionAccEtaEET->Fill(part->Energy(),part->Eta(),et);
		    fHistPionAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		    fHistPionAccEtaET->Fill(et,part->Eta());
		    fHistPionAccEtaE->Fill(part->Energy(),part->Eta());
		    fHistPionAccEtaPt->Fill(part->Pt(),part->Eta());													
		  }						
	      }
	    if (pdg->PdgCode() == fgKPlusCode || pdg->PdgCode() == fgKMinusCode)
	      {
		fKaontotET += et;
		fHistKaonEtaEET->Fill(part->Energy(),part->Eta(),et);
		fHistKaonEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		fHistKaonEtaE->Fill(part->Energy(),part->Eta());
		fHistKaonEtaET->Fill(et,part->Eta());
		fHistKaonEtaPt->Fill(part->Pt(),part->Eta());							
		// inside EMCal acceptance
		//if (IsInAcceptance(part,pdg,extParam)) 
		if (IsInAcceptance(aliPart)) 
		  {
		    fKaonAcctotET += et;
		    fHistKaonAccEtaEET->Fill(part->Energy(),part->Eta(),et);
		    fHistKaonAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		    fHistKaonAccEtaET->Fill(et,part->Eta());
		    fHistKaonAccEtaE->Fill(part->Energy(),part->Eta());
		    fHistKaonAccEtaPt->Fill(part->Pt(),part->Eta());													
		  }						
	      }
	    if (pdg->PdgCode() == fgMuPlusCode || pdg->PdgCode() == fgMuMinusCode)
	      {
		fMuontotET += et;
		fHistMuonEtaEET->Fill(part->Energy(),part->Eta(),et);
		fHistMuonEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		fHistMuonEtaET->Fill(et,part->Eta());
		fHistMuonEtaE->Fill(part->Energy(),part->Eta());
		fHistMuonEtaPt->Fill(part->Pt(),part->Eta());										
		// inside EMCal acceptance
		//if (IsInAcceptance(part,pdg,extParam)) 
		if (IsInAcceptance(aliPart)) 
		  {
		    fMuonAcctotET += et;
		    fHistMuonAccEtaEET->Fill(part->Energy(),part->Eta(),et);
		    fHistMuonAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		    fHistMuonAccEtaET->Fill(et,part->Eta());
		    fHistMuonAccEtaE->Fill(part->Energy(),part->Eta());
		    fHistMuonAccEtaPt->Fill(part->Pt(),part->Eta());																
		  }						
	      }
	    if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
	      {				
		if (stack->IsPhysicalPrimary(iPart))
		  {//Marcelo - isn't this redundant?  Isn't this inside an if statement already?
		    fElectrontotET += et;						
		    fHistElectronEtaEET->Fill(part->Energy(),part->Eta(),et);
		    fHistElectronEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		    fHistElectronEtaET->Fill(et,part->Eta());
		    fHistElectronEtaE->Fill(part->Energy(),part->Eta());
		    fHistElectronEtaPt->Fill(part->Pt(),part->Eta());	
						
		    // inside EMCal acceptance
		    //if (IsInAcceptance(part,pdg,extParam)) 
		    if (IsInAcceptance(aliPart)) 
		      {
			fElectronAcctotET += et;
			fHistElectronAccEtaEET->Fill(part->Energy(),part->Eta(),et);
			fHistElectronAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			fHistElectronAccEtaET->Fill(et,part->Eta());
			fHistElectronAccEtaE->Fill(part->Energy(),part->Eta());
			fHistElectronAccEtaPt->Fill(part->Pt(),part->Eta());							
		      }							
		  }
		else if (!fGeoUt->IsInEMCAL(part->Vx(),part->Vy(),part->Vz())) 
		  {//Marcelo - are we sure we know what this is doing?  How sensitive is this to geometry?
		    if (IsMotherPrimaryGamma(stack,iPartMom,pdgMom))
		      {
			fHistConvElectronEtaEET->Fill(part->Energy(),part->Eta(),et);
			fHistConvElectronEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			fHistConvElectronEtaET->Fill(et,part->Eta());
			fHistConvElectronEtaE->Fill(part->Energy(),part->Eta());
			fHistConvElectronEtaPt->Fill(part->Pt(),part->Eta());							
			fConvElectrontotET += et;
							
			// gamma mother is inside EMCal acceptance
			//if (IsInAcceptance(partMom,pdgMom)) 
			if (IsInAcceptance(aliPart)) 
			  {
			    fHistConvElectronAccEtaEET->Fill(part->Energy(),part->Eta(),et);
			    fHistConvElectronAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			    fHistConvElectronAccEtaET->Fill(et,part->Eta());
			    fHistConvElectronAccEtaE->Fill(part->Energy(),part->Eta());
			    fHistConvElectronAccEtaPt->Fill(part->Pt(),part->Eta());							
			    fConvElectronAcctotET += et;
			  }							
		      }		
		    else if (IsMotherPrimaryElectron(stack,iPartMom,pdgMom))
		      {
			fHistScatElectronEtaEET->Fill(part->Energy(),part->Eta(),et);
			fHistScatElectronEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			fHistScatElectronEtaET->Fill(et,part->Eta());
			fHistScatElectronEtaE->Fill(part->Energy(),part->Eta());
			fHistScatElectronEtaPt->Fill(part->Pt(),part->Eta());							
			fScatElectrontotET += et;

			// inside EMCal acceptance - does it work?
			//if (IsInAcceptance(part,pdg,extParam)) 
			if (IsInAcceptance(aliPart)) 
			  {
			    fHistScatElectronAccEtaEET->Fill(part->Energy(),part->Eta(),et);
			    fHistScatElectronAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			    fHistScatElectronAccEtaET->Fill(et,part->Eta());
			    fHistScatElectronAccEtaE->Fill(part->Energy(),part->Eta());
			    fHistScatElectronAccEtaPt->Fill(part->Pt(),part->Eta());							
			    fScatElectronAcctotET += et;
			  }							
		      }		
		  }

		// few checks
		if ((stack->IsPhysicalPrimary(iPart)) || (!fGeoUt->IsInEMCAL(part->Vx(),part->Vy(),part->Vz())))
		  {//Marcelo - ...isn't this redundant?
		    if (pdgMom)
		      fHistElectronFirstMother->Fill(pdgMom->PdgCode());
		    fHistElectronFirstMotherXY->Fill(part->Vx(),part->Vy());					
		    fHistElectronNDaughters->Fill(nPartDaughters);
						
		    iPartDaughter = part->GetLastDaughter();
		    if ((iPartDaughter>=0) && (iPartDaughter < nStackTracks))
		      {
			partDaughter = stack->Particle(iPartDaughter);
			if (partDaughter)
			  {
			    pdgDaugther = partDaughter->GetPDG(0);
			    if (pdgDaugther) {
			      fHistElectronDaughters->Fill(pdgDaugther->PdgCode());	
			      fHistElectronDaughtersXY->Fill(partDaughter->Vx(),partDaughter->Vy());
			    }
			  }
		      }
						
		    // inside EMCal acceptance
		    //if (IsInAcceptance(part,pdg,extParam)) 
		    if (IsInAcceptance(aliPart)) 
		      {//Marcelo - should all three of the lines below be in the if statement?
			if (pdgMom) fHistElectronFirstMotherAcc->Fill(pdgMom->PdgCode());
			fHistElectronFirstMotherXYAcc->Fill(part->Vx(),part->Vy());					
			fHistElectronNDaughtersAcc->Fill(nPartDaughters);
							
			iPartDaughter = part->GetLastDaughter();
			if ((iPartDaughter>=0) && (iPartDaughter < nStackTracks))
			  {
			    partDaughter = stack->Particle(iPartDaughter);
			    if (partDaughter)
			      {
				pdgDaugther = partDaughter->GetPDG(0);
				if (pdgDaugther) {
				  fHistElectronDaughtersAcc->Fill(pdgDaugther->PdgCode());	
				  fHistElectronDaughtersXYAcc->Fill(partDaughter->Vx(),partDaughter->Vy());
				}
			      }
			  }
		      }
		  }
	      } // end of if electron

	    // some neutrals also
	    if (pdg->PdgCode() == fgNeutronCode || pdg->PdgCode() == fgAntiNeutronCode)
	      {
		fHistNeutronEtaEET->Fill(part->Energy(),part->Eta(),et);
		fHistNeutronEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		fHistNeutronEtaET->Fill(et,part->Eta());							
		fHistNeutronEtaE->Fill(part->Energy(),part->Eta());
		fHistNeutronEtaPt->Fill(part->Pt(),part->Eta());		
		fNeutrontotET += et;
					
		// inside EMCal acceptance
		//if (IsInAcceptance(part,pdg)) 
		if (IsInAcceptance(aliPart)) 
		  {
		    fHistNeutronAccEtaEET->Fill(part->Energy(),part->Eta(),et);
		    fHistNeutronAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		    fHistNeutronAccEtaET->Fill(et,part->Eta());
		    fHistNeutronAccEtaE->Fill(part->Energy(),part->Eta());
		    fHistNeutronAccEtaPt->Fill(part->Pt(),part->Eta());													
		    fNeutronAcctotET += et;
		  }						
		if(pdg->PdgCode() == fgNeutronCode)
		  {
		    fNeutronEt += et;
		  }
		if(pdg->PdgCode() == fgAntiNeutronCode)
		  {
		    fAntiNeutronEt += et;
		  }
	      }
				
	    if(pdg->PdgCode() == fgGammaCode)
	      {
		if (stack->IsPhysicalPrimary(iPart))
		  {
		    fHistGammaEtaEET->Fill(part->Energy(),part->Eta(),et);
		    fHistGammaEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		    fHistGammaEtaET->Fill(et,part->Eta());
		    fHistGammaEtaE->Fill(part->Energy(),part->Eta());
		    fHistGammaEtaPt->Fill(part->Pt(),part->Eta());	
		    fGammatotET += et;
						
		    if (IsGammaConversion(stack, part, pdg))
		      {										
			fHistConvGammaEtaEET->Fill(part->Energy(),part->Eta(),et);
			fHistConvGammaEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			fHistConvGammaEtaET->Fill(et,part->Eta());
			fHistConvGammaEtaE->Fill(part->Energy(),part->Eta());
			fHistConvGammaEtaPt->Fill(part->Pt(),part->Eta());									
			fConvGammatotET += et;
		      }
		    else
		      {										
			fHistNonConvGammaEtaEET->Fill(part->Energy(),part->Eta(),et);
			fHistNonConvGammaEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			fHistNonConvGammaEtaET->Fill(et,part->Eta());
			fHistNonConvGammaEtaE->Fill(part->Energy(),part->Eta());
			fHistNonConvGammaEtaPt->Fill(part->Pt(),part->Eta());									
			fNonConvGammatotET += et;
		      }
						
		    Bool_t inAcc=kFALSE;
		    // inside EMCal acceptance
		    //if (IsInAcceptance(part,pdg)) 
		    if (IsInAcceptance(aliPart)) 
		      {
			//Printf("phi(1) = %f, eta(1) = %f",part->Phi(),part->Eta());
			inAcc = kTRUE;
							
			fHistGammaAccEtaEET->Fill(part->Energy(),part->Eta(),et);
			fHistGammaAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			fHistGammaAccEtaET->Fill(et,part->Eta());
			fHistGammaAccEtaE->Fill(part->Energy(),part->Eta());
			fHistGammaAccEtaPt->Fill(part->Pt(),part->Eta());		
			fGammaAcctotET += et;
		      }
						
		    if (IsInAcceptance(part,pdg)) 
		      {
			if (IsGammaConversion(stack, part, pdg))
			  {	
			    if (inAcc)
			      Printf("phi(1) = %f, eta(1) = %f",part->Phi(),part->Eta());

			    fHistConvGammaAccEtaEET->Fill(part->Energy(),part->Eta(),et);
			    fHistConvGammaAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			    fHistConvGammaAccEtaET->Fill(et,part->Eta());
			    fHistConvGammaAccEtaE->Fill(part->Energy(),part->Eta());
			    fHistConvGammaAccEtaPt->Fill(part->Pt(),part->Eta());									
			    fConvGammaAcctotET += et;
			  }
			else
			  {			
			    /*
			      if (!inAcc)
			      {
			      Printf("phi(2) = %f, eta(2) = %f, pt = %f",TMath::RadToDeg()*part->Phi(),part->Eta(),part->Pt());
								
			      Int_t iPartDaughter = part->GetLastDaughter();
									
			      if ((iPartDaughter>=0) && (iPartDaughter < nStackTracks))
			      {
			      TParticle *partDaughter = stack->Particle(iPartDaughter);
			      if (partDaughter)
			      {
			      TParticlePDG *pdgDaugther = partDaughter->GetPDG(0);
			      if (pdgDaugther) 
			      {
			      Double_t decayR = sqrt(pow(partDaughter->Vx(),2)+pow(partDaughter->Vy(),2));
			      Printf("radius = %f, daughter pid = %d",decayR,pdgDaugther->PdgCode());
			      }
			      }
			      }
									
			      for (int i=0;i<aliPart->GetNumberOfTrackReferences();i++)
			      {
			      AliTrackReference* aliTrkRef = aliPart->GetTrackReference(i);
										
			      if (aliTrkRef)
			      {
			      Printf("det id = %d, x=%f, y=%f, z=%f", aliTrkRef->DetectorId(),aliTrkRef->X(),aliTrkRef->Y(),aliTrkRef->Z());
			      }
			      }
									
			      }
			    */
								
			    fHistNonConvGammaAccEtaEET->Fill(part->Energy(),part->Eta(),et);
			    fHistNonConvGammaAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			    fHistNonConvGammaAccEtaET->Fill(et,part->Eta());
			    fHistNonConvGammaAccEtaE->Fill(part->Energy(),part->Eta());
			    fHistNonConvGammaAccEtaPt->Fill(part->Pt(),part->Eta());									
			    fNonConvGammaAcctotET += et;
			  }
		      }											
		  }
		else if (!fGeoUt->IsInEMCAL(part->Vx(),part->Vy(),part->Vz())) 
		  {
		    if (IsMotherPrimaryElectron(stack,iPartMom,pdgMom))
		      {
			fHistAnnihGammaEtaEET->Fill(part->Energy(),part->Eta(),et);
			fHistAnnihGammaEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			fHistAnnihGammaEtaET->Fill(et,part->Eta());
			fHistAnnihGammaEtaE->Fill(part->Energy(),part->Eta());
			fHistAnnihGammaEtaPt->Fill(part->Pt(),part->Eta());		
			fAnnihGammatotET += et;
							
			// inside EMCal acceptance
			//if (IsInAcceptance(part,pdg)) 
			if (IsInAcceptance(aliPart)) 
			  {
			    fHistAnnihGammaAccEtaEET->Fill(part->Energy(),part->Eta(),et);
			    fHistAnnihGammaAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			    fHistAnnihGammaAccEtaET->Fill(et,part->Eta());
			    fHistAnnihGammaAccEtaE->Fill(part->Energy(),part->Eta());
			    fHistAnnihGammaAccEtaPt->Fill(part->Pt(),part->Eta());												
			    fAnnihGammaAcctotET += et;
			  }
		      }
		    else if (IsMotherPrimaryGamma(stack,iPartMom,pdgMom))
		      {
			fHistScatGammaEtaEET->Fill(part->Energy(),part->Eta(),et);
			fHistScatGammaEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			fHistScatGammaEtaET->Fill(et,part->Eta());
			fHistScatGammaEtaE->Fill(part->Energy(),part->Eta());
			fHistScatGammaEtaPt->Fill(part->Pt(),part->Eta());		
			fScatGammatotET += et;
							
			// inside EMCal acceptance
			//if (IsInAcceptance(part,pdg)) 
			if (IsInAcceptance(aliPart)) 
			  {
			    fHistScatGammaAccEtaEET->Fill(part->Energy(),part->Eta(),et);
			    fHistScatGammaAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			    fHistScatGammaAccEtaET->Fill(et,part->Eta());
			    fHistScatGammaAccEtaE->Fill(part->Energy(),part->Eta());
			    fHistScatGammaAccEtaPt->Fill(part->Pt(),part->Eta());												
			    fScatGammaAcctotET += et;
			  }
		      }
		  }

		// few checks
		if ((stack->IsPhysicalPrimary(iPart)) || (!fGeoUt->IsInEMCAL(part->Vx(),part->Vy(),part->Vz())))
		  {
		    if (pdgMom)
		      fHistGammaFirstMother->Fill(pdgMom->PdgCode());
		    fHistGammaFirstMotherXY->Fill(part->Vx(),part->Vy());					
		    fHistGammaNDaughters->Fill(nPartDaughters);
						
		    iPartDaughter = part->GetLastDaughter();
		    if ((iPartDaughter>=0) && (iPartDaughter < nStackTracks))
		      {
			partDaughter = stack->Particle(iPartDaughter);
			if (partDaughter)
			  {
			    pdgDaugther = partDaughter->GetPDG(0);
			    if (pdgDaugther) {
			      fHistGammaDaughters->Fill(pdgDaugther->PdgCode());	
			      fHistGammaDaughtersXY->Fill(partDaughter->Vx(),partDaughter->Vy());
									
			      if (stack->IsPhysicalPrimary(iPart))
				{
				  if (IsGammaConversion(stack, part, pdg))
				    {										
				      fHistConvGammaDaughtersXY->Fill(partDaughter->Vx(),partDaughter->Vy());
				    }
				  else
				    {
				      fHistNonConvGammaDaughtersXY->Fill(partDaughter->Vx(),partDaughter->Vy());
				    }
				}
			    }
			  }
		      }
						
		    // inside EMCal acceptance
		    //if (IsInAcceptance(part,pdg)) 
		    if (IsInAcceptance(aliPart)) 
		      {
			if (pdgMom)//Marcelo - again, should this be all three?
			  fHistGammaFirstMotherAcc->Fill(pdgMom->PdgCode());
			fHistGammaFirstMotherXYAcc->Fill(part->Vx(),part->Vy());					
			fHistGammaNDaughtersAcc->Fill(nPartDaughters);
							
			iPartDaughter = part->GetLastDaughter();
			if ((iPartDaughter>=0) && (iPartDaughter < nStackTracks))
			  {
			    partDaughter = stack->Particle(iPartDaughter);
			    if (partDaughter)
			      {
				pdgDaugther = partDaughter->GetPDG(0);
				if (pdgDaugther) {
				  fHistGammaDaughtersAcc->Fill(pdgDaugther->PdgCode());	
				  fHistGammaDaughtersXYAcc->Fill(partDaughter->Vx(),partDaughter->Vy());
										
				  if (stack->IsPhysicalPrimary(iPart))
				    {
				      if (IsGammaConversion(stack, part, pdg))
					{										
					  fHistConvGammaDaughtersXYAcc->Fill(partDaughter->Vx(),partDaughter->Vy());
					}
				      else
					{
					  fHistNonConvGammaDaughtersXYAcc->Fill(partDaughter->Vx(),partDaughter->Vy());
					}
				    }
										
				}
			      }
			  }						
		      }
		  }
	      } // end of if gamma
				 	
	    // Neutral particles
	    if (TMath::Abs(pdg->Charge() - fCuts->GetMonteCarloNeutralParticle()) <1e-3 )
	      {
		//fNeutralMultiplicity++;
		fTotNeutralEt += et;
					
		// inside EMCal acceptance
		//if (IsInAcceptance(part,pdg)) 
		if (IsInAcceptance(aliPart)) 
		  {
		    fTotNeutralEtAcc += et;
		    //fTotEtAcc += et;										
		  }
	      } // end of neutral particles block
	    //Charged particles
	    else if (TMath::Abs( pdg->Charge() - fCuts->GetMonteCarloNeutralParticle())>1e-3 )
	      {
		//fChargedMultiplicity++;
		fTotChargedEt += et;
					
		// inside EMCal acceptance
		//if (IsInAcceptance(part,pdg,extParam)) 
		if (IsInAcceptance(aliPart)) 
		  {
		    fTotChargedEtAcc += et;
		    //fTotEtAcc += et;
		  } // inside EMCal acceptance
					
		//if (TrackHitsCalo(extParam)) // magnetic field info not filled?
		//{
		//	if (pdg->Charge() > 0) fHistPhivsPtPos->Fill(part->Phi(),part->Pt());
		//	else if (pdg->Charge() < 0) fHistPhivsPtNeg->Fill(part->Phi(), part->Pt());
		//}
	      } // end of charged particles block
	  } // end of is primary
	else // not a primary
	  {
	    if (pdgMom)
	      {
		if (pdgMom->PdgCode() == fgK0SCode)
		  {
		    fHistK0EtaEET->Fill(part->Energy(),part->Eta(),et);
		    fHistK0EtaPtET->Fill(part->Pt(),part->Eta(),et);							
		    fHistK0EtaET->Fill(et,part->Eta());							
		    fHistK0EtaE->Fill(part->Energy(),part->Eta());
		    fHistK0EtaPt->Fill(part->Pt(),part->Eta());		
		    fK0totET += et;
		  }
					
		if (pdgMom->PdgCode() == fgLambdaCode || pdgMom->PdgCode() == fgAntiLambdaCode)
		  {
		    fHistLambdaEtaEET->Fill(part->Energy(),part->Eta(),et);
		    fHistLambdaEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		    fHistLambdaEtaET->Fill(et,part->Eta());							
		    fHistLambdaEtaE->Fill(part->Energy(),part->Eta());
		    fHistLambdaEtaPt->Fill(part->Pt(),part->Eta());		
		    fLambdatotET += et;
		  }
	      }
								
	    if (!fGeoUt->IsInEMCAL(part->Vx(),part->Vy(),part->Vz())) // exclude secondaries from interactions inside the EMCal
	      {
		if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
		  {
		    fHistNPPElectronEtaEET->Fill(part->Energy(),part->Eta(),et);
		    fHistNPPElectronEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		    fHistNPPElectronEtaET->Fill(et,part->Eta());
		    fHistNPPElectronEtaE->Fill(part->Energy(),part->Eta());
		    fHistNPPElectronEtaPt->Fill(part->Pt(),part->Eta());	
		    fNPPElectrontotET += et;
						
		    // inside EMCal acceptance
		    //if (IsInAcceptance(part,pdg,extParam)) 
		    if (IsInAcceptance(aliPart)) 
		      {
			fHistNPPElectronAccEtaEET->Fill(part->Energy(),part->Eta(),et);
			fHistNPPElectronAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			fHistNPPElectronAccEtaE->Fill(part->Energy(),part->Eta());
			fHistNPPElectronAccEtaPt->Fill(part->Pt(),part->Eta());							
		      }							
						
		    // few checks
		    fHistNPPElectronFirstMother->Fill(pdgMom->PdgCode());
		    fHistNPPElectronFirstMotherXY->Fill(part->Vx(),part->Vy());					
		    fHistNPPElectronNDaughters->Fill(nPartDaughters);
						
		    iPartDaughter = part->GetLastDaughter();
		    if ((iPartDaughter>=0) && (iPartDaughter < nStackTracks))
		      {
			partDaughter = stack->Particle(iPartDaughter);
			if (partDaughter)
			  {
			    pdgDaugther = partDaughter->GetPDG(0);
			    if (pdgDaugther) {
			      fHistNPPElectronDaughters->Fill(pdgDaugther->PdgCode());	
			      fHistNPPElectronDaughtersXY->Fill(partDaughter->Vx(),partDaughter->Vy());
			    }
			  }
		      }
						
		    // inside EMCal acceptance
		    //if (IsInAcceptance(part,pdg,extParam)) 
		    if (IsInAcceptance(aliPart)) 
		      {
			fHistNPPElectronFirstMotherAcc->Fill(pdgMom->PdgCode());
			fHistNPPElectronFirstMotherXYAcc->Fill(part->Vx(),part->Vy());					
			fHistNPPElectronNDaughtersAcc->Fill(nPartDaughters);
						
			iPartDaughter = part->GetLastDaughter();
			if ((iPartDaughter>=0) && (iPartDaughter < nStackTracks))
			  {
			    partDaughter = stack->Particle(iPartDaughter);
			    if (partDaughter)
			      {
				pdgDaugther = partDaughter->GetPDG(0);
				if (pdgDaugther) {
				  fHistNPPElectronDaughtersAcc->Fill(pdgDaugther->PdgCode());	
				  fHistNPPElectronDaughtersXYAcc->Fill(partDaughter->Vx(),partDaughter->Vy());
				}
			      }
			  }
		      }						
						
		  } // end of if electron
					
		if(pdg->PdgCode() == fgGammaCode)
		  {
		    fHistNPPGammaEtaEET->Fill(part->Energy(),part->Eta(),et);
		    fHistNPPGammaEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		    fHistNPPGammaEtaET->Fill(et,part->Eta());
		    fHistNPPGammaEtaE->Fill(part->Energy(),part->Eta());
		    fHistNPPGammaEtaPt->Fill(part->Pt(),part->Eta());			
		    fNPPGammatotET += et;
						
		    if (pdgMom)
		      {	
			if (pdgMom->PdgCode() == fgPi0Code)
			  {
			    fHistNPPPi0GammaEtaEET->Fill(part->Energy(),part->Eta(),et);
			    fHistNPPPi0GammaEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			    fHistNPPPi0GammaEtaET->Fill(et,part->Eta());
			    fHistNPPPi0GammaEtaE->Fill(part->Energy(),part->Eta());
			    fHistNPPPi0GammaEtaPt->Fill(part->Pt(),part->Eta());									
			    fNPPPi0GammatotET += et;
			  }
		      }
						
		    // inside EMCal acceptance
		    //if (IsInAcceptance(part,pdg)) 
		    if (IsInAcceptance(aliPart)) 
		      {
			fHistNPPGammaAccEtaEET->Fill(part->Energy(),part->Eta(),et);
			fHistNPPGammaAccEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			fHistNPPGammaAccEtaE->Fill(part->Energy(),part->Eta());
			fHistNPPGammaAccEtaPt->Fill(part->Pt(),part->Eta());																		
		      }											
						
		    // few checks
		    if (pdgMom)
		      fHistNPPGammaFirstMother->Fill(pdgMom->PdgCode());
		    fHistNPPGammaFirstMotherXY->Fill(part->Vx(),part->Vy());					
		    fHistNPPGammaNDaughters->Fill(nPartDaughters);
						
		    iPartDaughter = part->GetLastDaughter();
		    if ((iPartDaughter>=0)  && (iPartDaughter < nStackTracks))
		      {
			partDaughter = stack->Particle(iPartDaughter);
			if (partDaughter)
			  {
			    pdgDaugther = partDaughter->GetPDG(0);
			    if (pdgDaugther) {
			      fHistNPPGammaDaughters->Fill(pdgDaugther->PdgCode());	
			      fHistNPPGammaDaughtersXY->Fill(partDaughter->Vx(),partDaughter->Vy());
			    }
			  }
		      }	
						
		    // inside EMCal acceptance
		    //if (IsInAcceptance(part,pdg)) 
		    if (IsInAcceptance(aliPart)) 
		      {
			fHistNPPGammaFirstMotherAcc->Fill(pdgMom->PdgCode());
			fHistNPPGammaFirstMotherXYAcc->Fill(part->Vx(),part->Vy());					
			fHistNPPGammaNDaughtersAcc->Fill(nPartDaughters);
						
			iPartDaughter = part->GetLastDaughter();
			if ((iPartDaughter>=0) && (iPartDaughter < nStackTracks))
			  {
			    partDaughter = stack->Particle(iPartDaughter);
			    if (partDaughter)
			      {
				pdgDaugther = partDaughter->GetPDG(0);
				if (pdgDaugther) {
				  fHistNPPGammaDaughtersAcc->Fill(pdgDaugther->PdgCode());	
				  fHistNPPGammaDaughtersXYAcc->Fill(partDaughter->Vx(),partDaughter->Vy());
				}
			      }
			  }
		      }
						
		  } // end of gamma
	      }
	  } // end of NOT a primary
      } // end of eta cut (Inside ALICE central barrel acceptance)
		
    if (extParam)
      delete extParam;
		
  }// end of loop over TParticles
  fTotEt = fTotChargedEt + fTotNeutralEt;
  fTotEtAcc = fTotChargedEtAcc + fTotNeutralEtAcc;   
	
  fTotElectrontotET = fElectrontotET + fConvElectrontotET + fScatElectrontotET;
  fTotElectronAcctotET = fElectronAcctotET + fConvElectronAcctotET + fScatElectronAcctotET;
  fTotGammatotET = fGammatotET + fAnnihGammatotET + fScatGammatotET;
  fTotGammaAcctotET = fGammaAcctotET + fAnnihGammaAcctotET + fScatGammaAcctotET;
  fTotEMtotET = fTotElectrontotET + fTotGammatotET;
  fTotEMAcctotET = fTotElectronAcctotET + fTotGammaAcctotET;
  fTotNPPEMtotET = fNPPElectrontotET + fNPPGammatotET;
  fTotChargedtotET = fMuontotET + fPiontotET + fKaontotET + fProtontotET;
  fTotChargedAcctotET = fMuonAcctotET + fPionAcctotET + fKaonAcctotET + fProtonAcctotET;
  fTotNeutraltotET = fNeutrontotET + fK0totET + fLambdatotET;
  fTotaltotET = fTotEMtotET + fTotNPPEMtotET + fTotChargedtotET + fTotNeutraltotET;
  fTotalAcctotET = fTotEMAcctotET + fTotChargedAcctotET;
	
  //FillHistograms();
	
  fHistPrimtotET->Fill(fPrimtotET);
  fHistPrimAcctotET->Fill(fPrimAcctotET);
	
  fHistElectrontotET->Fill(fElectrontotET);
  fHistElectronAcctotET->Fill(fElectronAcctotET);
  fHistConvElectrontotET->Fill(fConvElectrontotET);
  fHistConvElectronAcctotET->Fill(fConvElectronAcctotET);
  fHistScatElectrontotET->Fill(fScatElectrontotET);
  fHistScatElectronAcctotET->Fill(fScatElectronAcctotET);
	
  fHistTotElectrontotET->Fill(fTotElectrontotET);
  fHistTotElectronAcctotET->Fill(fTotElectronAcctotET);
	
  fHistGammatotET->Fill(fGammatotET);
  fHistGammaAcctotET->Fill(fGammaAcctotET);
  fHistAnnihGammatotET->Fill(fAnnihGammatotET);
  fHistAnnihGammaAcctotET->Fill(fAnnihGammaAcctotET);
  fHistScatGammatotET->Fill(fScatGammatotET);
  fHistScatGammaAcctotET->Fill(fScatGammaAcctotET);
	
  fHistTotGammatotET->Fill(fTotGammatotET);
  fHistTotGammaAcctotET->Fill(fTotGammaAcctotET);
	
  fHistTotEMtotET->Fill(fTotEMtotET);
  fHistTotEMAcctotET->Fill(fTotEMAcctotET);
	
  fHistConvGammatotET->Fill(fConvGammatotET);
  fHistNonConvGammatotET->Fill(fNonConvGammatotET);
  fHistConvGammaAcctotET->Fill(fConvGammaAcctotET);
  fHistNonConvGammaAcctotET->Fill(fNonConvGammaAcctotET);

  fHistNPPElectrontotET->Fill(fNPPElectrontotET);
  fHistNPPGammatotET->Fill(fNPPGammatotET);
	
  fHistTotNPPEMtotET->Fill(fTotNPPEMtotET);

  fHistNPPPi0GammatotET->Fill(fNPPPi0GammatotET);

  fHistMuontotET->Fill(fMuontotET); 
  fHistMuonAcctotET->Fill(fMuonAcctotET); 
  fHistPiontotET->Fill(fPiontotET); 
  fHistPionAcctotET->Fill(fPionAcctotET); 
  fHistKaontotET->Fill(fKaontotET); 
  fHistKaonAcctotET->Fill(fKaonAcctotET); 
  fHistProtontotET->Fill(fProtontotET); 
  fHistProtonAcctotET->Fill(fProtonAcctotET); 

  fHistTotChargedtotET->Fill(fTotChargedtotET);
  fHistTotChargedAcctotET->Fill(fTotChargedAcctotET);
	
  fHistNeutrontotET->Fill(fNeutrontotET); 
  fHistNeutronAcctotET->Fill(fNeutronAcctotET); 
  fHistK0totET->Fill(fK0totET); 
  fHistLambdatotET->Fill(fNeutrontotET); 

  fHistTotNeutraltotET->Fill(fTotNeutraltotET);
	
  fHistTotaltotET->Fill(fTotaltotET);
  fHistTotalAcctotET->Fill(fTotalAcctotET);
	 
  return 0;    
}

Int_t AliAnalysisEmEtMonteCarlo::AnalyseEvent(AliVEvent* ev,AliVEvent* ev2)
{ // analyse MC and real event info
  if(!ev || !ev2){//Marcelo - should use AliError
    Printf("ERROR: Event does not exist");   
    return 0;
  }
	
  AliMCEvent *mcEvent = dynamic_cast<AliMCEvent*>(ev);
  AliESDEvent *realEvent = dynamic_cast<AliESDEvent*>(ev2);

  if(!fGeoUt){
    fGeoUt = AliEMCALGeometry::GetInstance("EMCAL_FIRSTYEARV1");//new AliEMCALGeometry("EMCAL_FIRSTYEAR","EMCAL");
    AliInfo("Creating new AliEMCALGeometry");
  }
  //fGeoUt = new AliEMCALGeometry("EMCAL_FIRSTYEAR","EMCAL");
  //fGeoUt->SetMisalMatrix(realEvent->GetEMCALMatrix(0),0);

  ResetEventValues();
  AnalyseEvent(ev);
	
  AliStack *stack = mcEvent->Stack();
  if (!stack)
    {
      Printf("ERROR: Could not get stack");
      return 0;
    }	
	
  Int_t nStackTracks = stack->GetNtrack();
	
  // get all emcal clusters
  TRefArray* caloClusters = new TRefArray();
  realEvent->GetEMCALClusters( caloClusters );
	
  Int_t nCluster = caloClusters->GetEntries();
	
  Float_t pos[3] = {0};
  TVector3 caloPos(0,0,0);
  TVector3 trackPos(0,0,0);
	
  // loop the clusters
  for (int iCluster = 0; iCluster < nCluster; iCluster++ ) 
    {
      AliESDCaloCluster* caloCluster = ( AliESDCaloCluster* )caloClusters->At( iCluster );
      Float_t caloE = caloCluster->E();
      caloCluster->GetPosition(pos);		
      caloPos.SetXYZ(pos[0],pos[1],pos[2]);
		
      UInt_t iPart = (UInt_t)TMath::Abs(caloCluster->GetLabel());
      TParticle *part  = stack->Particle(iPart);
		
      if (!part)
        {//Marcelo -- use AliError
	  Printf("No MC particle %d", iCluster);
	  continue;
        }
		
      TParticlePDG *pdg = part->GetPDG(0);
		
      TParticle *partMom = 0;
      TParticlePDG *pdgMom = 0;
		
      Int_t nPartDaughters = part->GetNDaughters();
      TParticle *partDaughter = 0;
      TParticlePDG *pdgDaugther = 0;		
		
      if (!pdg)
        {//Marcelo -- use AliError
	  Printf("ERROR-2: Could not get particle PDG %d", iPart);
	  continue;
        }		
		
      Int_t iPartMom = part->GetMother(0);
      Int_t iPartDaughter = 0;
		
      if ((iPartMom>=0) && (iPartMom < nStackTracks))
	{
	  partMom = stack->Particle(iPartMom);
	  pdgMom = partMom->GetPDG(0);
	}			
		
      // find the track associated to this MC particle
      TObjArray* list = fEsdtrackCutsITSTPC->GetAcceptedTracks(realEvent);
      Int_t nGoodTracks = list->GetEntries();
      Bool_t trackFound = kFALSE;
      Bool_t trackProjected = kFALSE;
      Float_t res = 0;
      AliESDtrack *track = 0;
      AliEMCALTrack *emcTrack = 0;
      AliExternalTrackParam* extParamTPart = 0;
      AliESDtrack *esdTPart = 0;
      AliEMCALTrack *emcTPart = 0;
			
      // find corresponding track
      for (Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++)
	{
	  track = dynamic_cast<AliESDtrack*> (list->At(iTrack));
	  if (!track)
	    {//Marcelo -use AliError
	      Printf("ERROR: Could not get track %d", iTrack);
	      continue;
	    }
	  else
	    {
	      UInt_t label = (UInt_t)TMath::Abs(track->GetLabel());
	      if (label == iPart)
		{
		  trackFound = kTRUE;
		  emcTrack = new AliEMCALTrack(*track);
		  fHistAllPtRecPtMC->Fill(part->Pt(),track->Pt());
					
		  if (GetTrackProjection(emcTrack,trackPos,caloPos)) 
		    {
		      trackProjected = kTRUE;
		      res = sqrt(pow(trackPos.Phi()-caloPos.Phi(),2)+pow(trackPos.Eta()-caloPos.Eta(),2));
		    }
		  else
		    res = -1.;
					
		  break;					
		}
	    }
	}			

      if (!trackFound)
	{
	  track = 0;
	  emcTrack = 0;
	  res = -2.;
	}
			
      //create an external track param for projection
      extParamTPart = CreateExternalTrackParam(part);
		
      // create esd and emcal tracks out of TParticle (used for projection)
      esdTPart = new AliESDtrack(part);
		
      if (esdTPart && extParamTPart) 
	{
	  esdTPart->SetOuterParam(extParamTPart,0);
	  emcTPart = new AliEMCALTrack(*esdTPart);
	}
		
      // few checks
      // compare MC and Rec energies for all particles
      fHistAllERecEMC->Fill(part->Energy(),caloE);
      //Marcelo - doesn't it make sense to change this so that we use the tracks matched by the official code?		
      if (TMath::Abs( pdg->Charge() - fCuts->GetMonteCarloNeutralParticle()) > 1e-3)
	{
	  //Printf("calo.Phi = %f, calo.Eta = %f \n", caloPos.Phi(), caloPos.Eta());
			
	  if (trackProjected)
	    {
	      //Printf("good track.Phi = %f, track.Eta = %f  \n", trackPos.Phi(), trackPos.Eta());
	      fHistChargedRes->Fill(trackPos.Phi()-caloPos.Phi(),trackPos.Eta()-caloPos.Eta());
	    }
			
	  if (GetTrackProjection(emcTPart,trackPos,caloPos)) 
	    {
	      fHistChargedRes2->Fill(trackPos.Phi()-caloPos.Phi(),trackPos.Eta()-caloPos.Eta());
	      //Printf("track.Phi = %f, track.Eta = %f  \n", trackPos.Phi(), trackPos.Eta());
	    }
			
	  if (GetTrackProjection(extParamTPart,trackPos)) 
	    {
	      fHistChargedRes3->Fill(trackPos.Phi()-caloPos.Phi(),trackPos.Eta()-caloPos.Eta());
	      //Printf("track.Phi = %f, track.Eta = %f  \n", trackPos.Phi(), trackPos.Eta());
	    }
			
	} 
      else if (TMath::Abs(pdg->Charge() - fCuts->GetMonteCarloNeutralParticle()) < 1e-3 )
	{
	  fHistNeutralRes->Fill(part->Phi()-caloPos.Phi(),part->Eta()-caloPos.Eta());			
	}
		
      if(pdg->PdgCode() == fgGammaCode)
	{
	  // compare MC and Rec energies for gammas
	  fHistGammaERecEMC->Fill(part->Energy(),caloE);		
	  fHistGammaRes->Fill(part->Phi()-caloPos.Phi(),part->Eta()-caloPos.Eta());
	}			
		
      if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
	{
	  // compare MC and Rec energies for electrons
	  fHistElectronERecEMC->Fill(part->Energy(),caloE);
	  if (GetTrackProjection(extParamTPart,trackPos))
	    { 
	      fHistElectronRes->Fill(trackPos.Phi()-caloPos.Phi(),trackPos.Eta()-caloPos.Eta());
	    }			
	}
		
      // calculate ET
      Double_t et = CalcET(part,pdg);
      Double_t etDep = CalcETDep(caloE,part,pdg);
		
      // Check if it is a primary particle
      if (IsPrimary(stack,iPart,pdg,iPartMom,pdgMom))
	{				
	  if (TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloSingleChargedParticle())<1e-3 && TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloNeutralParticle())<1e-3) continue;
				
	  if (stack->IsPhysicalPrimary(iPart))
	    {
	      fHistPrimRecEtaEET->Fill(part->Energy(),part->Eta(),et);
	      fHistPrimRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
	      fHistPrimRecEtaET->Fill(et,part->Eta());
	      fPrimRectotET += et;

	      fHistPrimRecEtaEDepETDep->Fill(part->Energy(),part->Eta(),etDep);
	      fHistPrimRecEtaPtETDep->Fill(part->Pt(),part->Eta(),etDep);							
	      fHistPrimRecEtaETDep->Fill(etDep,part->Eta());
	      fPrimRectotETDep += etDep;
	    }
			 			
	  if(pdg->PdgCode() == fgGammaCode)
	    {			
	      if (stack->IsPhysicalPrimary(iPart))
		{
		  fHistGammaRecEtaEET->Fill(part->Energy(),part->Eta(),et);
		  fHistGammaRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		  fHistGammaRecEtaET->Fill(et,part->Eta());
		  fHistGammaRecEtaE->Fill(part->Energy(),part->Eta());
		  fHistGammaRecEtaPt->Fill(part->Pt(),part->Eta());	
		  fGammaRectotET += et;
		}
	      else if (!fGeoUt->IsInEMCAL(part->Vx(),part->Vy(),part->Vz())) 
		{
		  if (IsMotherPrimaryElectron(stack,iPartMom,pdgMom))
		    {
		      fHistAnnihGammaRecEtaEET->Fill(part->Energy(),part->Eta(),et);
		      fHistAnnihGammaRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		      fHistAnnihGammaRecEtaET->Fill(et,part->Eta());
		      fHistAnnihGammaRecEtaE->Fill(part->Energy(),part->Eta());
		      fHistAnnihGammaRecEtaPt->Fill(part->Pt(),part->Eta());		
		      fAnnihGammaRectotET += et;
		    }					
		  else if (IsMotherPrimaryGamma(stack,iPartMom,pdgMom))
		    {
		      fHistScatGammaRecEtaEET->Fill(part->Energy(),part->Eta(),et);
		      fHistScatGammaRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		      fHistScatGammaRecEtaET->Fill(et,part->Eta());
		      fHistScatGammaRecEtaE->Fill(part->Energy(),part->Eta());
		      fHistScatGammaRecEtaPt->Fill(part->Pt(),part->Eta());		
		      fScatGammaRectotET += et;
		    }					
		}
				
	      // few checks
	      if (pdgMom)
		fHistGammaFirstMotherRec->Fill(pdgMom->PdgCode());
	      fHistGammaFirstMotherXYRec->Fill(part->Vx(),part->Vy());					
	      fHistGammaNDaughtersRec->Fill(nPartDaughters);
				
	      iPartDaughter = part->GetLastDaughter();
	      if ((iPartDaughter>=0)  && (iPartDaughter < nStackTracks))
		{
		  partDaughter = stack->Particle(iPartDaughter);
		  if (partDaughter)
		    {
		      pdgDaugther = partDaughter->GetPDG(0);
		      if (pdgDaugther) {
			fHistGammaDaughtersRec->Fill(pdgDaugther->PdgCode());	
			fHistGammaDaughtersXYRec->Fill(partDaughter->Vx(),partDaughter->Vy());
							
			if (stack->IsPhysicalPrimary(iPart))
			  {
			    if (IsGammaConversion(stack, part, pdg))
			      {										
				fHistConvGammaDaughtersXYRec->Fill(partDaughter->Vx(),partDaughter->Vy());
			      }
			    else
			      {
				fHistNonConvGammaDaughtersXYRec->Fill(partDaughter->Vx(),partDaughter->Vy());
			      }
			  }							
		      }
		    }
		}
	    } // gamma
			
	  if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
	    {
	      if (stack->IsPhysicalPrimary(iPart))
		{
		  fHistElectronRecEtaEET->Fill(part->Energy(),part->Eta(),et);
		  fHistElectronRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		  fHistElectronRecEtaET->Fill(et,part->Eta());
		  fHistElectronRecEtaE->Fill(part->Energy(),part->Eta());
		  fHistElectronRecEtaPt->Fill(part->Pt(),part->Eta());	
		  fElectronRectotET += et;
		}
	      else if (!fGeoUt->IsInEMCAL(part->Vx(),part->Vy(),part->Vz())) 
		{
		  if (IsMotherPrimaryGamma(stack,iPartMom,pdgMom))
		    {
		      fHistConvElectronRecEtaEET->Fill(part->Energy(),part->Eta(),et);
		      fHistConvElectronRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		      fHistConvElectronRecEtaET->Fill(et,part->Eta());
		      fHistConvElectronRecEtaE->Fill(part->Energy(),part->Eta());
		      fHistConvElectronRecEtaPt->Fill(part->Pt(),part->Eta());							
		      fConvElectronRectotET += et;
		    }		
		  else if (IsMotherPrimaryElectron(stack,iPartMom,pdgMom))
		    {
		      fHistScatElectronRecEtaEET->Fill(part->Energy(),part->Eta(),et);
		      fHistScatElectronRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		      fHistScatElectronRecEtaET->Fill(et,part->Eta());
		      fHistScatElectronRecEtaE->Fill(part->Energy(),part->Eta());
		      fHistScatElectronRecEtaPt->Fill(part->Pt(),part->Eta());							
		      fScatElectronRectotET += et;
		    }		
		}
				
	      // few checks
	      if (pdgMom)
		fHistElectronFirstMotherRec->Fill(pdgMom->PdgCode());
	      fHistElectronFirstMotherXYRec->Fill(part->Vx(),part->Vy());					
	      fHistElectronNDaughtersRec->Fill(nPartDaughters);
				
	      iPartDaughter = part->GetLastDaughter();
	      if ((iPartDaughter>=0) && (iPartDaughter < nStackTracks))
		{
		  partDaughter = stack->Particle(iPartDaughter);
		  if (partDaughter)
		    {
		      pdgDaugther = partDaughter->GetPDG(0);
		      if (pdgDaugther) {
			fHistElectronDaughtersRec->Fill(pdgDaugther->PdgCode());	
			fHistElectronDaughtersXYRec->Fill(partDaughter->Vx(),partDaughter->Vy());
		      }
		    }
		}
	    } // electrons
			
	  if (pdg->PdgCode() == fgMuPlusCode || pdg->PdgCode() == fgMuMinusCode)
	    {
	      fHistMuonRecEtaEET->Fill(part->Energy(),part->Eta(),et);
	      fHistMuonRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
	      fHistMuonRecEtaET->Fill(et,part->Eta());
	      fHistMuonRecEtaE->Fill(part->Energy(),part->Eta());
	      fHistMuonRecEtaPt->Fill(part->Pt(),part->Eta());	
	      fMuonRectotET += et;
				
	      fHistMuonRecEtaEDepETDep->Fill(caloE,part->Eta(),etDep);
	      fHistMuonRecEtaPtETDep->Fill(part->Pt(),part->Eta(),etDep);							
	      fHistMuonRecEtaETDep->Fill(etDep,part->Eta());
	      fMuonRectotETDep += etDep;

	      if (trackProjected)
		{
		  fHistMuonRecResEET->Fill(part->Energy(),res,et);
		  fHistMuonRecResPtET->Fill(part->Pt(),res,et);							
		  fHistMuonRecResE->Fill(part->Energy(),res);
		  fHistMuonRecResPt->Fill(part->Pt(),res);	
		  fHistMuonRecResEDepETDep->Fill(caloE,res,etDep);
		  fHistMuonRecResPtETDep->Fill(part->Pt(),res,etDep);	
					
		  if ((res>0.) && (res<fResCut))
		    {
		      fHistMuonMatchEtaEET->Fill(part->Energy(),part->Eta(),et);
		      fHistMuonMatchEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		      fHistMuonMatchEtaET->Fill(et,part->Eta());
		      fHistMuonMatchEtaE->Fill(part->Energy(),part->Eta());
		      fHistMuonMatchEtaPt->Fill(part->Pt(),part->Eta());	
		      fMuonMatchtotET += et;
						
		      fHistMuonMatchEtaEDepETDep->Fill(caloE,part->Eta(),etDep);
		      fHistMuonMatchEtaPtETDep->Fill(part->Pt(),part->Eta(),etDep);							
		      fHistMuonMatchEtaETDep->Fill(etDep,part->Eta());
		      fMuonMatchtotETDep += etDep;						
		    }
		}
	    }
			
	  if (pdg->PdgCode() == fgPiPlusCode || pdg->PdgCode() == fgPiMinusCode)
	    {
	      fHistPionRecEtaEET->Fill(part->Energy(),part->Eta(),et);
	      fHistPionRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
	      fHistPionRecEtaET->Fill(et,part->Eta());
	      fHistPionRecEtaE->Fill(part->Energy(),part->Eta());
	      fHistPionRecEtaPt->Fill(part->Pt(),part->Eta());	
	      fPionRectotET += et;

	      fHistPionRecEtaEDepETDep->Fill(caloE,part->Eta(),etDep);
	      fHistPionRecEtaPtETDep->Fill(part->Pt(),part->Eta(),etDep);							
	      fHistPionRecEtaETDep->Fill(etDep,part->Eta());
	      fPionRectotETDep += etDep;

	      if (trackProjected)
		{
		  fHistPionRecResEET->Fill(part->Energy(),res,et);
		  fHistPionRecResPtET->Fill(part->Pt(),res,et);							
		  fHistPionRecResE->Fill(part->Energy(),res);
		  fHistPionRecResPt->Fill(part->Pt(),res);	
		  fHistPionRecResEDepETDep->Fill(caloE,res,etDep);
		  fHistPionRecResPtETDep->Fill(part->Pt(),res,etDep);							

		  if ((res>0.) && (res<fResCut))
		    {
		      fHistPionMatchEtaEET->Fill(part->Energy(),part->Eta(),et);
		      fHistPionMatchEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		      fHistPionMatchEtaET->Fill(et,part->Eta());
		      fHistPionMatchEtaE->Fill(part->Energy(),part->Eta());
		      fHistPionMatchEtaPt->Fill(part->Pt(),part->Eta());	
		      fPionMatchtotET += et;
						
		      fHistPionMatchEtaEDepETDep->Fill(caloE,part->Eta(),etDep);
		      fHistPionMatchEtaPtETDep->Fill(part->Pt(),part->Eta(),etDep);							
		      fHistPionMatchEtaETDep->Fill(etDep,part->Eta());
		      fPionMatchtotETDep += etDep;
		    }
		}
	    }			
			
	  if (pdg->PdgCode() == fgKPlusCode || pdg->PdgCode() == fgKMinusCode)
	    {
	      fHistKaonRecEtaEET->Fill(part->Energy(),part->Eta(),et);
	      fHistKaonRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
	      fHistKaonRecEtaET->Fill(et,part->Eta());
	      fHistKaonRecEtaE->Fill(part->Energy(),part->Eta());
	      fHistKaonRecEtaPt->Fill(part->Pt(),part->Eta());	
	      fKaonRectotET += et;

	      fHistKaonRecEtaEDepETDep->Fill(caloE,part->Eta(),etDep);
	      fHistKaonRecEtaPtETDep->Fill(part->Pt(),part->Eta(),etDep);							
	      fHistKaonRecEtaETDep->Fill(etDep,part->Eta());
	      fKaonRectotETDep += etDep;

	      if (trackProjected)
		{
		  fHistKaonRecResEET->Fill(part->Energy(),res,et);
		  fHistKaonRecResPtET->Fill(part->Pt(),res,et);							
		  fHistKaonRecResE->Fill(part->Energy(),res);
		  fHistKaonRecResPt->Fill(part->Pt(),res);	

		  fHistKaonRecResEDepETDep->Fill(caloE,res,etDep);
		  fHistKaonRecResPtETDep->Fill(part->Pt(),res,etDep);							

		  if ((res>0.) && (res<fResCut))
		    {
		      fHistKaonMatchEtaEET->Fill(part->Energy(),part->Eta(),et);
		      fHistKaonMatchEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		      fHistKaonMatchEtaET->Fill(et,part->Eta());
		      fHistKaonMatchEtaE->Fill(part->Energy(),part->Eta());
		      fHistKaonMatchEtaPt->Fill(part->Pt(),part->Eta());	
		      fKaonMatchtotET += et;
						
		      fHistKaonMatchEtaEDepETDep->Fill(caloE,part->Eta(),etDep);
		      fHistKaonMatchEtaPtETDep->Fill(part->Pt(),part->Eta(),etDep);							
		      fHistKaonMatchEtaETDep->Fill(etDep,part->Eta());
		      fKaonMatchtotETDep += etDep;
		    }
		}
	    }
	
	  if (pdg->PdgCode() == fgProtonCode || pdg->PdgCode() == fgAntiProtonCode)
	    {
	      fHistProtonRecEtaEET->Fill(part->Energy(),part->Eta(),et);
	      fHistProtonRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
	      fHistProtonRecEtaET->Fill(et,part->Eta());
	      fHistProtonRecEtaE->Fill(part->Energy(),part->Eta());
	      fHistProtonRecEtaPt->Fill(part->Pt(),part->Eta());	
	      fProtonRectotET += et;

	      fHistProtonRecEtaEDepETDep->Fill(caloE,part->Eta(),etDep);
	      fHistProtonRecEtaPtETDep->Fill(part->Pt(),part->Eta(),etDep);							
	      fHistProtonRecEtaETDep->Fill(etDep,part->Eta());
	      fProtonRectotETDep += etDep;

	      if (trackProjected)
		{
		  fHistProtonRecResEET->Fill(part->Energy(),res,et);
		  fHistProtonRecResPtET->Fill(part->Pt(),res,et);							
		  fHistProtonRecResE->Fill(part->Energy(),res);
		  fHistProtonRecResPt->Fill(part->Pt(),res);	

		  fHistProtonRecResEDepETDep->Fill(caloE,res,etDep);
		  fHistProtonRecResPtETDep->Fill(part->Pt(),res,etDep);							

		  if ((res>0.) && (res<fResCut))
		    {
		      fHistProtonMatchEtaEET->Fill(part->Energy(),part->Eta(),et);
		      fHistProtonMatchEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		      fHistProtonMatchEtaET->Fill(et,part->Eta());
		      fHistProtonMatchEtaE->Fill(part->Energy(),part->Eta());
		      fHistProtonMatchEtaPt->Fill(part->Pt(),part->Eta());	
		      fProtonMatchtotET += et;
						
		      fHistProtonMatchEtaEDepETDep->Fill(caloE,part->Eta(),etDep);
		      fHistProtonMatchEtaPtETDep->Fill(part->Pt(),part->Eta(),etDep);							
		      fHistProtonMatchEtaETDep->Fill(etDep,part->Eta());
		      fProtonMatchtotETDep += etDep;
		    }
		}				
	    }
			
	  if (pdg->PdgCode() == fgNeutronCode || pdg->PdgCode() == fgAntiNeutronCode)
	    {
	      fHistNeutronRecEtaEET->Fill(part->Energy(),part->Eta(),et);
	      fHistNeutronRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
	      fHistNeutronRecEtaET->Fill(et,part->Eta());
	      fHistNeutronRecEtaE->Fill(part->Energy(),part->Eta());
	      fHistNeutronRecEtaPt->Fill(part->Pt(),part->Eta());	
	      fNeutronRectotET += et;

	      fHistNeutronRecEtaEDepETDep->Fill(caloE,part->Eta(),etDep);
	      fHistNeutronRecEtaPtETDep->Fill(part->Pt(),part->Eta(),etDep);							
	      fHistNeutronRecEtaETDep->Fill(etDep,part->Eta());
	      fNeutronRectotETDep += etDep;
	    }
			
	  if (emcTrack)
	    delete emcTrack;
	  if (esdTPart)
	    delete esdTPart;
	  if (emcTPart)
	    delete emcTPart;
	  if (extParamTPart)
	    delete extParamTPart;
	} // end of primary tracks
      else // not a primary
	{
	  if (pdgMom)
	    {
	      if (pdgMom->PdgCode() == fgK0SCode)
		{
		  fHistK0RecEtaEET->Fill(part->Energy(),part->Eta(),et);
		  fHistK0RecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		  fHistK0RecEtaET->Fill(et,part->Eta());
		  fHistK0RecEtaE->Fill(part->Energy(),part->Eta());
		  fHistK0RecEtaPt->Fill(part->Pt(),part->Eta());	
		  fK0RectotET += et;
					
		  fHistK0RecEtaEDepETDep->Fill(caloE,part->Eta(),etDep);
		  fHistK0RecEtaPtETDep->Fill(part->Pt(),part->Eta(),etDep);							
		  fHistK0RecEtaETDep->Fill(etDep,part->Eta());					
		  fK0RectotETDep += etDep;					
		}
				
	      if (pdgMom->PdgCode() == fgLambdaCode || pdgMom->PdgCode() == fgAntiLambdaCode)
		{
		  fHistLambdaRecEtaEET->Fill(part->Energy(),part->Eta(),et);
		  fHistLambdaRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		  fHistLambdaRecEtaET->Fill(et,part->Eta());
		  fHistLambdaRecEtaE->Fill(part->Energy(),part->Eta());
		  fHistLambdaRecEtaPt->Fill(part->Pt(),part->Eta());	
		  fLambdaRectotET += et;

		  fHistLambdaRecEtaEDepETDep->Fill(caloE,part->Eta(),etDep);
		  fHistLambdaRecEtaPtETDep->Fill(part->Pt(),part->Eta(),etDep);							
		  fHistLambdaRecEtaETDep->Fill(etDep,part->Eta());					
		  fLambdaRectotETDep += etDep;
		}
	    }
			
	  if (!fGeoUt->IsInEMCAL(part->Vx(),part->Vy(),part->Vz())) // exclude secondaries from interactions inside the EMCal
	    {
	      if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
		{
		  fHistNPPElectronRecEtaEET->Fill(part->Energy(),part->Eta(),et);
		  fHistNPPElectronRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		  fHistNPPElectronRecEtaET->Fill(et,part->Eta());
		  fHistNPPElectronRecEtaE->Fill(part->Energy(),part->Eta());
		  fHistNPPElectronRecEtaPt->Fill(part->Pt(),part->Eta());	
		  fNPPElectronRectotET += et;
					
		  // few checks
		  fHistNPPElectronFirstMotherRec->Fill(pdgMom->PdgCode());
		  fHistNPPElectronFirstMotherXYRec->Fill(part->Vx(),part->Vy());					
		  fHistNPPElectronNDaughtersRec->Fill(nPartDaughters);
					
		  iPartDaughter = part->GetLastDaughter();
		  if ((iPartDaughter>=0) && (iPartDaughter < nStackTracks))
		    {
		      partDaughter = stack->Particle(iPartDaughter);
		      if (partDaughter)
			{
			  pdgDaugther = partDaughter->GetPDG(0);
			  if (pdgDaugther) {
			    fHistNPPElectronDaughtersRec->Fill(pdgDaugther->PdgCode());	
			    fHistNPPElectronDaughtersXYRec->Fill(partDaughter->Vx(),partDaughter->Vy());
			  }
			}
		    }
		} // end of if electron
				
	      if(pdg->PdgCode() == fgGammaCode)
		{
		  fHistNPPGammaRecEtaEET->Fill(part->Energy(),part->Eta(),et);
		  fHistNPPGammaRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
		  fHistNPPGammaRecEtaET->Fill(et,part->Eta());
		  fHistNPPGammaRecEtaE->Fill(part->Energy(),part->Eta());
		  fHistNPPGammaRecEtaPt->Fill(part->Pt(),part->Eta());			
		  fNPPGammaRectotET += et;
					
		  if (pdgMom)
		    {	
		      if (pdgMom->PdgCode() == fgPi0Code)
			{
			  fHistNPPPi0GammaRecEtaEET->Fill(part->Energy(),part->Eta(),et);
			  fHistNPPPi0GammaRecEtaPtET->Fill(part->Pt(),part->Eta(),et);							
			  fHistNPPPi0GammaRecEtaET->Fill(et,part->Eta());
			  fHistNPPPi0GammaRecEtaE->Fill(part->Energy(),part->Eta());
			  fHistNPPPi0GammaRecEtaPt->Fill(part->Pt(),part->Eta());									
			  fNPPPi0GammaRectotET += et;
			}
		    }
					
		  // few checks
		  if (pdgMom)
		    fHistNPPGammaFirstMotherRec->Fill(pdgMom->PdgCode());
		  fHistNPPGammaFirstMotherXYRec->Fill(part->Vx(),part->Vy());					
		  fHistNPPGammaNDaughtersRec->Fill(nPartDaughters);
					
		  iPartDaughter = part->GetLastDaughter();
		  if ((iPartDaughter>=0)  && (iPartDaughter < nStackTracks))
		    {
		      partDaughter = stack->Particle(iPartDaughter);
		      if (partDaughter)
			{
			  pdgDaugther = partDaughter->GetPDG(0);
			  if (pdgDaugther) {
			    fHistNPPGammaDaughtersRec->Fill(pdgDaugther->PdgCode());	
			    fHistNPPGammaDaughtersXYRec->Fill(partDaughter->Vx(),partDaughter->Vy());
			  }
			}
		    }	
		} // end of gamma
	    }
	} // end of NOT a primary
    } // end of loop over clusters	

  fTotElectronRectotET = fElectronRectotET + fConvElectronRectotET + fScatElectronRectotET;
  fTotGammaRectotET = fGammaRectotET + fAnnihGammaRectotET + fScatElectronRectotET;
  fTotEMRectotET = fTotElectronRectotET + fTotGammaRectotET;
  fTotNPPEMRectotET = fNPPElectronRectotET + fNPPGammaRectotET;
  fTotChargedRectotET = fMuonRectotET + fPionRectotET + fKaonRectotET + fProtonRectotET;
  fTotChargedRectotETDep = fMuonRectotETDep + fPionRectotETDep + fKaonRectotETDep + fProtonRectotETDep;
  fTotChargedMatchtotET = fMuonMatchtotET + fPionMatchtotET + fKaonMatchtotET + fProtonMatchtotET;
  fTotChargedMatchtotETDep = fMuonMatchtotETDep + fPionMatchtotETDep + fKaonMatchtotETDep + fProtonMatchtotETDep;
  fTotNeutralRectotET = fNeutronRectotET + fK0RectotET + fLambdaRectotET;
  fTotNeutralRectotETDep = fNeutronRectotETDep + fK0RectotETDep + fLambdaRectotETDep;
  fTotalRectotET = fTotEMRectotET + fTotNPPEMRectotET + fTotChargedRectotET + fTotNeutralRectotET;
  fTotalRectotETDep = fTotEMRectotET + fTotNPPEMRectotET + fTotChargedRectotETDep + fTotNeutralRectotETDep;
	
  fHistPrimRectotET->Fill(fPrimRectotET);
  fHistPrimRectotETDep->Fill(fPrimRectotETDep);

  fHistElectronRectotET->Fill(fElectronRectotET);
  fHistConvElectronRectotET->Fill(fConvElectronRectotET);
  fHistScatElectronRectotET->Fill(fScatElectronRectotET);
  fHistTotElectronRectotET->Fill(fTotElectronRectotET);
	
  fHistGammaRectotET->Fill(fGammaRectotET);
  fHistAnnihGammaRectotET->Fill(fAnnihGammaRectotET);
  fHistScatGammaRectotET->Fill(fScatGammaRectotET);
  fHistTotGammaRectotET->Fill(fTotGammaRectotET);
	
  fHistTotEMRectotET->Fill(fTotEMRectotET);
	
  fHistNPPElectronRectotET->Fill(fNPPElectronRectotET);
  fHistNPPGammaRectotET->Fill(fNPPGammaRectotET);
  fHistTotNPPEMRectotET->Fill(fTotNPPEMRectotET);
	
  fHistNPPPi0GammaRectotET->Fill(fNPPPi0GammaRectotET);

  fHistMuonRectotET->Fill(fMuonRectotET); 
  fHistMuonRectotETDep->Fill(fMuonRectotETDep); 
  fHistMuonMatchtotET->Fill(fMuonMatchtotET); 
  fHistMuonMatchtotETDep->Fill(fMuonMatchtotETDep); 
  fHistPionRectotET->Fill(fPionRectotET); 
  fHistPionRectotETDep->Fill(fPionRectotETDep); 
  fHistPionMatchtotET->Fill(fPionMatchtotET); 
  fHistPionMatchtotETDep->Fill(fPionMatchtotETDep); 
  fHistKaonRectotET->Fill(fKaonRectotET); 
  fHistKaonRectotETDep->Fill(fKaonRectotETDep); 
  fHistKaonMatchtotET->Fill(fKaonMatchtotET); 
  fHistKaonMatchtotETDep->Fill(fKaonMatchtotETDep); 
  fHistProtonRectotET->Fill(fProtonRectotET); 
  fHistProtonRectotETDep->Fill(fProtonRectotETDep); 
  fHistProtonMatchtotET->Fill(fProtonMatchtotET); 
  fHistProtonMatchtotETDep->Fill(fProtonMatchtotETDep); 
  fHistTotChargedRectotET->Fill(fTotChargedRectotET);
  fHistTotChargedRectotETDep->Fill(fTotChargedRectotETDep);
  fHistTotChargedMatchtotET->Fill(fTotChargedMatchtotET);
  fHistTotChargedMatchtotETDep->Fill(fTotChargedMatchtotETDep);
	
  fHistNeutronRectotET->Fill(fNeutronRectotET);
  fHistNeutronRectotETDep->Fill(fNeutronRectotETDep);
  fHistK0RectotET->Fill(fK0RectotET);
  fHistK0RectotETDep->Fill(fK0RectotETDep);
  fHistLambdaRectotET->Fill(fLambdaRectotET);
  fHistLambdaRectotETDep->Fill(fLambdaRectotETDep);
  fHistTotNeutralRectotET->Fill(fTotNeutralRectotET);
  fHistTotNeutralRectotETDep->Fill(fTotNeutralRectotETDep);
	
  fHistTotalRectotET->Fill(fTotalRectotET);
  fHistTotalRectotETDep->Fill(fTotalRectotETDep);
	
  //delete fGeoUt;
  delete caloClusters;//Marcelo - Christine - make as object & don't create new one each event
	
  return 0;    
}	

void AliAnalysisEmEtMonteCarlo::Init()
{ // init
  AliAnalysisEt::Init();
	
  fDetectorRadius = fCuts->GetGeometryEmcalDetectorRadius();
  fEtaCutAcc = fCuts->GetGeometryEmcalEtaAccCut();
  fPhiCutAccMax = fCuts->GetGeometryEmcalPhiAccMaxCut() * TMath::Pi()/180.;
  fPhiCutAccMin = fCuts->GetGeometryEmcalPhiAccMinCut() * TMath::Pi()/180.;
  fClusterEnergyCut = fCuts->GetReconstructedEmcalClusterEnergyCut();
  fSingleCellEnergyCut = fCuts->GetReconstructedEmcalSingleCellEnergyCut();
	
  fDetector = fCuts->GetDetectorEmcal();	
}

void AliAnalysisEmEtMonteCarlo::ResetEventValues()
{ // reset event values
  AliAnalysisEt::ResetEventValues();
	
  fPrimtotET = 0; fPrimAcctotET = 0; fPrimRectotET = 0; fPrimRectotETDep = 0;

  fElectrontotET = 0; fElectronAcctotET = 0; fElectronRectotET = 0;
  fConvElectrontotET = 0; fConvElectronAcctotET = 0; fConvElectronRectotET = 0; fScatElectrontotET = 0; fScatElectronAcctotET = 0; fScatElectronRectotET = 0;
  fTotElectrontotET = 0, fTotElectronAcctotET = 0, fTotElectronRectotET = 0;
	
  fGammatotET = 0; fGammaAcctotET = 0; fGammaRectotET = 0;
  fAnnihGammatotET = 0; fAnnihGammaAcctotET = 0; fAnnihGammaRectotET = 0; fScatGammatotET = 0; fScatGammaAcctotET = 0; fScatGammaRectotET = 0;
  fTotGammatotET = 0, fTotGammaAcctotET = 0, fTotGammaRectotET = 0;

  fTotEMtotET = 0, fTotEMAcctotET = 0, fTotEMRectotET = 0;

  fConvGammatotET = 0; fNonConvGammatotET = 0; fConvGammaAcctotET = 0; fNonConvGammaAcctotET = 0; fNPPPi0GammatotET = 0; fNPPPi0GammaRectotET = 0;
	
  fNPPElectrontotET = 0; fNPPElectronRectotET = 0; fNPPGammatotET = 0; fNPPGammaRectotET = 0;
  fTotNPPEMtotET = 0, fTotNPPEMRectotET = 0;

  fMuontotET = 0; fPiontotET = 0; fKaontotET = 0; fProtontotET = 0;
  fMuonAcctotET = 0; fPionAcctotET = 0; fKaonAcctotET = 0; fProtonAcctotET = 0;
  fMuonRectotET = 0; fMuonRectotETDep = 0; fPionRectotET = 0; fPionRectotETDep = 0; fKaonRectotET = 0; fKaonRectotETDep = 0; fProtonRectotET = 0; fProtonRectotETDep = 0;
  fMuonMatchtotET = 0; fMuonMatchtotETDep = 0; fPionMatchtotET = 0; fPionMatchtotETDep = 0; fKaonMatchtotET = 0; fKaonMatchtotETDep = 0; fProtonMatchtotET = 0; fProtonMatchtotETDep = 0;
  fTotChargedtotET = 0, fTotChargedAcctotET = 0, fTotChargedRectotET = 0, fTotChargedRectotETDep = 0, fTotChargedMatchtotET = 0, fTotChargedMatchtotETDep = 0;

  fNeutrontotET = 0; fNeutronAcctotET = 0; fNeutronRectotET = 0; fNeutronRectotETDep = 0;
  fK0totET = 0; fK0RectotET = 0; fK0RectotETDep = 0; fLambdatotET = 0; fLambdaRectotET = 0; fLambdaRectotETDep = 0;
  fTotNeutraltotET = 0, fTotNeutralRectotET = 0, fTotNeutralRectotETDep = 0;
	
  fTotaltotET = 0, fTotalAcctotET = 0, fTotalRectotET = 0, fTotalRectotETDep = 0;
}

void AliAnalysisEmEtMonteCarlo::CreateHistograms()
{ // histogram related additions
  //AliAnalysisEt::CreateHistograms();
	
  fHistPrimEtaEET = CreateEtaEHisto2D("fHistPrimEtaEET_","MC E_{T}, primary particles","E_{T}(GeV)");
  fHistPrimEtaPtET = CreateEtaPtHisto2D("fHistPrimEtaPtET_","MC E_{T}, primary particles","E_{T}(GeV)");
  fHistPrimEtaET = CreateEtaEtHisto2D("fHistPrimEtaET_","MC primary particles","#");
  TString histname = "fHistPrimtotET_" + fHistogramNameSuffix;
  fHistPrimtotET = new TH1F(histname.Data(),"total ET, primary particles",fgNumOfEBins, fgEAxis);
	
  fHistPrimAccEtaEET = CreateEtaEHisto2D("fHistPrimAccEtaEET_","MC E_{T}, primary particles","E_{T}(GeV)");
  fHistPrimAccEtaPtET = CreateEtaPtHisto2D("fHistPrimAccEtaPtET_","MC E_{T}, primary particles","E_{T}(GeV)");
  fHistPrimAccEtaET = CreateEtaEtHisto2D("fHistPrimAccEtaET_","MC primary particles","#");
  histname = "fHistPrimAcctotET_" + fHistogramNameSuffix;
  fHistPrimAcctotET = new TH1F(histname.Data(),"total ET, primary particles",fgNumOfEBins, fgEAxis);
	
  fHistPrimRecEtaEET = CreateEtaEHisto2D("fHistPrimRecEtaEET_","MC E_{T}, primary particles","E_{T}(GeV)");
  fHistPrimRecEtaPtET = CreateEtaPtHisto2D("fHistPrimRecEtaPtET_","MC E_{T}, primary particles","E_{T}(GeV)");
  fHistPrimRecEtaET = CreateEtaEtHisto2D("fHistPrimRecEtaET_","MC primary particles","#");
  histname = "fHistPrimRectotET_" + fHistogramNameSuffix;
  fHistPrimRectotET = new TH1F(histname.Data(),"total ET, primary particles",fgNumOfEBins, fgEAxis);

  fHistPrimRecEtaEDepETDep = CreateEtaEHisto2D("fHistPrimRecEtaEDepETDep_","MC E_{T}, primary particles","E_{T}(GeV)");
  fHistPrimRecEtaPtETDep = CreateEtaPtHisto2D("fHistPrimRecEtaPtETDep_","MC E_{T}, primary particles","E_{T}(GeV)");
  fHistPrimRecEtaETDep = CreateEtaEtHisto2D("fHistPrimRecEtaETDep_","MC primary particles","#");
  histname = "fHistPrimRectotETDep_" + fHistogramNameSuffix;
  fHistPrimRectotETDep = new TH1F(histname.Data(),"total ET, primary particles",fgNumOfEBins, fgEAxis);
	
  fHistElectronEtaEET = CreateEtaEHisto2D("fHistElectronEtaEET_","MC E_{T}, primary electrons","E_{T}(GeV)");
  fHistElectronEtaPtET = CreateEtaPtHisto2D("fHistElectronEtaPtET_","MC E_{T}, primary electrons","E_{T}(GeV)");
  fHistElectronEtaET = CreateEtaEtHisto2D("fHistElectronEtaET_","MC primary electrons","#");
  fHistElectronEtaE = CreateEtaEHisto2D("fHistElectronEtaE_","MC primary electrons","#");
  fHistElectronEtaPt = CreateEtaPtHisto2D("fHistElectronEtaPt_","MC primary electrons","#");
  histname = "fHistElectrontotET_" + fHistogramNameSuffix;
  fHistElectrontotET = new TH1F(histname.Data(),"total ET, MC primary electrons",fgNumOfEBins, fgEAxis);

  fHistConvElectronEtaEET = CreateEtaEHisto2D("fHistConvElectronEtaEET_","MC E_{T}, electrons from conversion","E_{T}(GeV)");
  fHistConvElectronEtaPtET = CreateEtaPtHisto2D("fHistConvElectronEtaPtET_","MC E_{T}, electrons from conversion","E_{T}(GeV)");
  fHistConvElectronEtaET = CreateEtaEtHisto2D("fHistConvElectronEtaET_","MC electrons from conversion","#");
  fHistConvElectronEtaE = CreateEtaEHisto2D("fHistConvElectronEtaE_","MC electrons from conversion","#");
  fHistConvElectronEtaPt = CreateEtaPtHisto2D("fHistConvElectronEtaPt_","MC electrons from conversion","#");
  histname = "fHistConvElectrontotET_" + fHistogramNameSuffix;
  fHistConvElectrontotET = new TH1F(histname.Data(),"total ET, MC electrons from conversion",fgNumOfEBins, fgEAxis);
	
  fHistScatElectronEtaEET = CreateEtaEHisto2D("fHistScatElectronEtaEET_","MC E_{T}, electrons from Scattering","E_{T}(GeV)");
  fHistScatElectronEtaPtET = CreateEtaPtHisto2D("fHistScatElectronEtaPtET_","MC E_{T}, electrons from Scattering","E_{T}(GeV)");
  fHistScatElectronEtaET = CreateEtaEtHisto2D("fHistScatElectronEtaET_","MC electrons from Scattering","#");
  fHistScatElectronEtaE = CreateEtaEHisto2D("fHistScatElectronEtaE_","MC electrons from Scattering","#");
  fHistScatElectronEtaPt = CreateEtaPtHisto2D("fHistScatElectronEtaPt_","MC electrons from Scattering","#");
  histname = "fHistScatElectrontotET_" + fHistogramNameSuffix;
  fHistScatElectrontotET = new TH1F(histname.Data(),"total ET, MC electrons from Scattering",fgNumOfEBins, fgEAxis);
	
  histname = "fHistTotElectrontotET_" + fHistogramNameSuffix;
  fHistTotElectrontotET = new TH1F(histname.Data(),"total ET, MC primary electrons",fgNumOfEBins, fgEAxis);

  fHistGammaEtaEET = CreateEtaEHisto2D("fHistGammaEtaEET_","MC E_{T}, primary gammas","E_{T}(GeV)"); 
  fHistGammaEtaPtET = CreateEtaPtHisto2D("fHistGammaEtaPtET_","MC E_{T}, primary gammas","E_{T}(GeV)"); 
  fHistGammaEtaET = CreateEtaEtHisto2D("fHistGammaEtaET_","MC primary gammas","#"); 
  fHistGammaEtaE = CreateEtaEHisto2D("fHistGammaEtaE_","MC primary gammas","#"); 
  fHistGammaEtaPt = CreateEtaPtHisto2D("fHistGammaEtaPt_","MC primary gammas","#"); 
  histname = "fHistGammatotET_" + fHistogramNameSuffix;
  fHistGammatotET = new TH1F(histname.Data(),"total ET, MC primary gammas",fgNumOfEBins, fgEAxis);

  fHistAnnihGammaEtaEET = CreateEtaEHisto2D("fHistAnnihGammaEtaEET_","MC E_{T}, Annihilation gammas","E_{T}(GeV)"); 
  fHistAnnihGammaEtaPtET = CreateEtaPtHisto2D("fHistAnnihGammaEtaPtET_","MC E_{T}, Annihilation gammas","E_{T}(GeV)"); 
  fHistAnnihGammaEtaET = CreateEtaEtHisto2D("fHistAnnihGammaEtaET_","MC Annihilation gammas","#"); 	
  fHistAnnihGammaEtaE = CreateEtaEHisto2D("fHistAnnihGammaEtaE_","MC Annihilation gammas","#"); 	
  fHistAnnihGammaEtaPt = CreateEtaPtHisto2D("fHistAnnihGammaEtaPt_","MC Annihilation gammas","#"); 
  histname = "fHistAnnihGammatotET_" + fHistogramNameSuffix;
  fHistAnnihGammatotET = new TH1F(histname.Data(),"total ET, MC Annihilation gammas",fgNumOfEBins, fgEAxis);
	
  fHistScatGammaEtaEET = CreateEtaEHisto2D("fHistScatGammaEtaEET_","MC E_{T}, Scattering gammas","E_{T}(GeV)"); 
  fHistScatGammaEtaPtET = CreateEtaPtHisto2D("fHistScatGammaEtaPtET_","MC E_{T}, Scattering gammas","E_{T}(GeV)"); 
  fHistScatGammaEtaET = CreateEtaEtHisto2D("fHistScatGammaEtaET_","MC Scattering gammas","#"); 	
  fHistScatGammaEtaE = CreateEtaEHisto2D("fHistScatGammaEtaE_","MC Scattering gammas","#"); 	
  fHistScatGammaEtaPt = CreateEtaPtHisto2D("fHistScatGammaEtaPt_","MC Scattering gammas","#"); 
  histname = "fHistScatGammatotET_" + fHistogramNameSuffix;
  fHistScatGammatotET = new TH1F(histname.Data(),"total ET, MC Scattering gammas",fgNumOfEBins, fgEAxis);
	
  fHistConvGammaEtaEET = CreateEtaEHisto2D("fHistConvGammaEtaEET_","MC E_{T}, non conversion primary gammas","E_{T}(GeV)"); 
  fHistConvGammaEtaPtET = CreateEtaPtHisto2D("fHistConvGammaEtaPtET_","MC E_{T}, non conversion primary gammas","E_{T}(GeV)"); 
  fHistConvGammaEtaET = CreateEtaEtHisto2D("fHistConvGammaEtaET_","MC non conversion primary gammas","#"); 
  fHistConvGammaEtaE = CreateEtaEHisto2D("fHistConvGammaEtaE_","MC non conversion primary gammas","#"); 
  fHistConvGammaEtaPt = CreateEtaPtHisto2D("fHistConvGammaEtaPt_","MC non conversion primary gammas","#"); 
  histname = "fHistConvGammatotET_" + fHistogramNameSuffix;
  fHistConvGammatotET = new TH1F(histname.Data(),"total ET, MC non conversion primary gammas",fgNumOfEBins, fgEAxis);
	
  fHistNonConvGammaEtaEET = CreateEtaEHisto2D("fHistNonConvGammaEtaEET_","MC E_{T}, non conversion primary gammas","E_{T}(GeV)"); 
  fHistNonConvGammaEtaPtET = CreateEtaPtHisto2D("fHistNonConvGammaEtaPtET_","MC E_{T}, non conversion primary gammas","E_{T}(GeV)"); 
  fHistNonConvGammaEtaET = CreateEtaEtHisto2D("fHistNonConvGammaEtaET_","MC non conversion primary gammas","#"); 
  fHistNonConvGammaEtaE = CreateEtaEHisto2D("fHistNonConvGammaEtaE_","MC non conversion primary gammas","#"); 
  fHistNonConvGammaEtaPt = CreateEtaPtHisto2D("fHistNonConvGammaEtaPt_","MC non conversion primary gammas","#"); 
  histname = "fHistNonConvGammatotET_" + fHistogramNameSuffix;
  fHistNonConvGammatotET = new TH1F(histname.Data(),"total ET, MC non conversion primary gammas",fgNumOfEBins, fgEAxis);
		
  histname = "fHistTotGammatotET_" + fHistogramNameSuffix;
  fHistTotGammatotET = new TH1F(histname.Data(),"total ET, MC primary gammas",fgNumOfEBins, fgEAxis);

  histname = "fHistTotEMtotET_" + fHistogramNameSuffix;
  fHistTotEMtotET = new TH1F(histname.Data(),"total electromagnetic ET",fgNumOfEBins, fgEAxis);

  fHistNPPElectronEtaEET = CreateEtaEHisto2D("fHistNPPElectronEtaEET_","MC E_{T}, non-primary electrons","E_{T}(GeV)");
  fHistNPPElectronEtaPtET = CreateEtaPtHisto2D("fHistNPPElectronEtaPtET_","MC E_{T}, non-primary electrons","E_{T}(GeV)");
  fHistNPPElectronEtaET = CreateEtaEtHisto2D("fHistNPPElectronEtaET_","MC non-primary electrons","#");
  fHistNPPElectronEtaE = CreateEtaEHisto2D("fHistNPPElectronEtaE_","MC non-primary electrons","#");
  fHistNPPElectronEtaPt = CreateEtaPtHisto2D("fHistNPPElectronEtaPt_","MC non-primary electrons","#");
  histname = "fHistNPPElectrontotET_" + fHistogramNameSuffix;
  fHistNPPElectrontotET = new TH1F(histname.Data(),"total ET, MC non-primary electrons",fgNumOfEBins, fgEAxis);
	
  fHistNPPGammaEtaEET = CreateEtaEHisto2D("fHistNPPGammaEtaEET_","MC E_{T}, non-primary gammas","E_{T}(GeV)"); 
  fHistNPPGammaEtaPtET = CreateEtaPtHisto2D("fHistNPPGammaEtaPtET_","MC E_{T}, non-primary gammas","E_{T}(GeV)"); 
  fHistNPPGammaEtaET = CreateEtaEtHisto2D("fHistNPPGammaEtaET_","MC non-primary gammas","#"); 
  fHistNPPGammaEtaE = CreateEtaEHisto2D("fHistNPPGammaEtaE_","MC non-primary gammas","#"); 
  fHistNPPGammaEtaPt = CreateEtaPtHisto2D("fHistNPPGammaEtaPt_","MC non-primary gammas","#"); 
  histname = "fHistNPPGammatotET_" + fHistogramNameSuffix;
  fHistNPPGammatotET = new TH1F(histname.Data(),"total ET, MC non-primary gammas",fgNumOfEBins, fgEAxis);
	
  histname = "fHistTotNPPEMtotET_" + fHistogramNameSuffix;
  fHistTotNPPEMtotET = new TH1F(histname.Data(),"total ET, MC non-primary electromagnetic",fgNumOfEBins, fgEAxis);
	
  fHistNPPPi0GammaEtaEET = CreateEtaEHisto2D("fHistNPPPi0GammaEtaEET_","MC E_{T}, non-primary gammas","E_{T}(GeV)"); 
  fHistNPPPi0GammaEtaPtET = CreateEtaPtHisto2D("fHistNPPPi0GammaEtaPtET_","MC E_{T}, non-primary gammas","E_{T}(GeV)"); 
  fHistNPPPi0GammaEtaET = CreateEtaEtHisto2D("fHistNPPPi0GammaEtaET_","MC non-primary gammas","#"); 
  fHistNPPPi0GammaEtaE = CreateEtaEHisto2D("fHistNPPPi0GammaEtaE_","MC non-primary gammas","#"); 
  fHistNPPPi0GammaEtaPt = CreateEtaPtHisto2D("fHistNPPPi0GammaEtaPt_","MC non-primary gammas","#"); 
  histname = "fHistNPPPi0GammatotET_" + fHistogramNameSuffix;
  fHistNPPPi0GammatotET = new TH1F(histname.Data(),"total ET, MC non-primary gammas",fgNumOfEBins, fgEAxis);
		
  fHistElectronAccEtaEET = CreateEtaEHisto2D("fHistElectronAccEtaEET_","MC E_{T}, primary electrons","E_{T}(GeV)");
  fHistElectronAccEtaPtET = CreateEtaPtHisto2D("fHistElectronAccEtaPtET_","MC E_{T}, primary electrons","E_{T}(GeV)");
  fHistElectronAccEtaET = CreateEtaEtHisto2D("fHistElectronAccEtaET_","MC primary electrons","#");
  fHistElectronAccEtaE = CreateEtaEHisto2D("fHistElectronAccEtaE_","MC primary electrons","#");
  fHistElectronAccEtaPt = CreateEtaPtHisto2D("fHistElectronAccEtaPt_","MC primary electrons","#");
  histname = "fHistElectronAcctotET_" + fHistogramNameSuffix;
  fHistElectronAcctotET = new TH1F(histname.Data(),"total ET, MC primary electrons",fgNumOfEBins, fgEAxis);
	
  fHistConvElectronAccEtaEET = CreateEtaEHisto2D("fHistConvElectronAccEtaEET_","MC E_{T}, electrons from conversion","E_{T}(GeV)");
  fHistConvElectronAccEtaPtET = CreateEtaPtHisto2D("fHistConvElectronAccEtaPtET_","MC E_{T}, electrons from conversion","E_{T}(GeV)");
  fHistConvElectronAccEtaET = CreateEtaEtHisto2D("fHistConvElectronAccEtaET_","MC electrons from conversion","#");
  fHistConvElectronAccEtaE = CreateEtaEHisto2D("fHistConvElectronAccEtaE_","MC electrons from conversion","#");
  fHistConvElectronAccEtaPt = CreateEtaPtHisto2D("fHistConvElectronAccEtaPt_","MC electrons from conversion","#");
  histname = "fHistConvElectronAcctotET_" + fHistogramNameSuffix;
  fHistConvElectronAcctotET = new TH1F(histname.Data(),"total ET, MC electrons from conversion",fgNumOfEBins, fgEAxis);
	
  fHistScatElectronAccEtaEET = CreateEtaEHisto2D("fHistScatElectronAccEtaEET_","MC E_{T}, electrons from Scattering","E_{T}(GeV)");
  fHistScatElectronAccEtaPtET = CreateEtaPtHisto2D("fHistScatElectronAccEtaPtET_","MC E_{T}, electrons from Scattering","E_{T}(GeV)");
  fHistScatElectronAccEtaET = CreateEtaEtHisto2D("fHistScatElectronAccEtaET_","MC electrons from Scattering","#");
  fHistScatElectronAccEtaE = CreateEtaEHisto2D("fHistScatElectronAccEtaE_","MC electrons from Scattering","#");
  fHistScatElectronAccEtaPt = CreateEtaPtHisto2D("fHistScatElectronAccEtaPt_","MC electrons from Scattering","#");
  histname = "fHistScatElectronAcctotET_" + fHistogramNameSuffix;
  fHistScatElectronAcctotET = new TH1F(histname.Data(),"total ET, MC electrons from Scattering",fgNumOfEBins, fgEAxis);

  histname = "fHistTotElectronAcctotET_" + fHistogramNameSuffix;
  fHistTotElectronAcctotET = new TH1F(histname.Data(),"total ET, MC primary electrons",fgNumOfEBins, fgEAxis);

  fHistGammaAccEtaEET = CreateEtaEHisto2D("fHistGammaAccEtaEET_","MC E_{T}, primary gammas","E_{T}(GeV)"); 
  fHistGammaAccEtaPtET = CreateEtaPtHisto2D("fHistGammaAccEtaPtET_","MC E_{T}, primary gammas","E_{T}(GeV)"); 
  fHistGammaAccEtaET = CreateEtaEtHisto2D("fHistGammaAccEtaET_","MC primary gammas","#"); 
  fHistGammaAccEtaE = CreateEtaEHisto2D("fHistGammaAccEtaE_","MC primary gammas","#"); 
  fHistGammaAccEtaPt = CreateEtaPtHisto2D("fHistGammaAccEtaPt_","MC primary gammas","#"); 
  histname = "fHistGammaAcctotET_" + fHistogramNameSuffix;
  fHistGammaAcctotET = new TH1F(histname.Data(),"total ET, MC primary gammas",fgNumOfEBins, fgEAxis);
	
  fHistAnnihGammaAccEtaEET = CreateEtaEHisto2D("fHistAnnihGammaAccEtaEET_","MC E_{T}, Annihilation gammas","E_{T}(GeV)"); 
  fHistAnnihGammaAccEtaPtET = CreateEtaPtHisto2D("fHistAnnihGammaAccEtaPtET_","MC E_{T}, Annihilation gammas","E_{T}(GeV)"); 
  fHistAnnihGammaAccEtaET = CreateEtaEtHisto2D("fHistAnnihGammaAccEtaET_","MC Annihilation gammas","#"); 	
  fHistAnnihGammaAccEtaE = CreateEtaEHisto2D("fHistAnnihGammaAccEtaE_","MC Annihilation gammas","#"); 	
  fHistAnnihGammaAccEtaPt = CreateEtaPtHisto2D("fHistAnnihGammaAccEtaPt_","MC Annihilation gammas","#"); 
  histname = "fHistAnnihGammaAcctotET_" + fHistogramNameSuffix;
  fHistAnnihGammaAcctotET = new TH1F(histname.Data(),"total ET, MC Annihilation gammas",fgNumOfEBins, fgEAxis);
	
  fHistScatGammaAccEtaEET = CreateEtaEHisto2D("fHistScatGammaAccEtaEET_","MC E_{T}, Scattering gammas","E_{T}(GeV)"); 
  fHistScatGammaAccEtaPtET = CreateEtaPtHisto2D("fHistScatGammaAccEtaPtET_","MC E_{T}, Scattering gammas","E_{T}(GeV)"); 
  fHistScatGammaAccEtaET = CreateEtaEtHisto2D("fHistScatGammaAccEtaET_","MC Scattering gammas","#"); 	
  fHistScatGammaAccEtaE = CreateEtaEHisto2D("fHistScatGammaAccEtaE_","MC Scattering gammas","#"); 	
  fHistScatGammaAccEtaPt = CreateEtaPtHisto2D("fHistScatGammaAccEtaPt_","MC Scattering gammas","#"); 
  histname = "fHistScatGammaAcctotET_" + fHistogramNameSuffix;
  fHistScatGammaAcctotET = new TH1F(histname.Data(),"total ET, MC Scattering gammas",fgNumOfEBins, fgEAxis);
	
  fHistConvGammaAccEtaEET = CreateEtaEHisto2D("fHistConvGammaAccEtaEET_","MC E_{T}, non conversion primary gammas","E_{T}(GeV)"); 
  fHistConvGammaAccEtaPtET = CreateEtaPtHisto2D("fHistConvGammaAccEtaPtET_","MC E_{T}, non conversion primary gammas","E_{T}(GeV)"); 
  fHistConvGammaAccEtaET = CreateEtaEtHisto2D("fHistConvGammaAccEtaET_","MC non conversion primary gammas","#"); 
  fHistConvGammaAccEtaE = CreateEtaEHisto2D("fHistConvGammaAccEtaE_","MC non conversion primary gammas","#"); 
  fHistConvGammaAccEtaPt = CreateEtaPtHisto2D("fHistConvGammaAccEtaPt_","MC non conversion primary gammas","#"); 
  histname = "fHistConvGammaAcctotET_" + fHistogramNameSuffix;
  fHistConvGammaAcctotET = new TH1F(histname.Data(),"total ET, MC non conversion primary gammas",fgNumOfEBins, fgEAxis);
	
  fHistNonConvGammaAccEtaEET = CreateEtaEHisto2D("fHistNonConvGammaAccEtaEET_","MC E_{T}, non conversion primary gammas","E_{T}(GeV)"); 
  fHistNonConvGammaAccEtaPtET = CreateEtaPtHisto2D("fHistNonConvGammaAccEtaPtET_","MC E_{T}, non conversion primary gammas","E_{T}(GeV)"); 
  fHistNonConvGammaAccEtaET = CreateEtaEtHisto2D("fHistNonConvGammaAccEtaET_","MC non conversion primary gammas","#"); 
  fHistNonConvGammaAccEtaE = CreateEtaEHisto2D("fHistNonConvGammaAccEtaE_","MC non conversion primary gammas","#"); 
  fHistNonConvGammaAccEtaPt = CreateEtaPtHisto2D("fHistNonConvGammaAccEtaPt_","MC non conversion primary gammas","#"); 
  histname = "fHistNonConvGammaAcctotET_" + fHistogramNameSuffix;
  fHistNonConvGammaAcctotET = new TH1F(histname.Data(),"total ET, MC non conversion primary gammas",fgNumOfEBins, fgEAxis);
	
  histname = "fHistTotGammaAcctotET_" + fHistogramNameSuffix;
  fHistTotGammaAcctotET = new TH1F(histname.Data(),"total ET, MC primary gammas",fgNumOfEBins, fgEAxis);
	
  histname = "fHistTotEMAcctotET_" + fHistogramNameSuffix;
  fHistTotEMAcctotET = new TH1F(histname.Data(),"total electromagnetic ET",fgNumOfEBins, fgEAxis);

  fHistNPPElectronAccEtaEET = CreateEtaEHisto2D("fHistNPPElectronAccEtaEET_","MC E_{T}, non-primary electrons","E_{T}(GeV)");
  fHistNPPElectronAccEtaPtET = CreateEtaPtHisto2D("fHistNPPElectronAccEtaPtET_","MC E_{T}, non-primary electrons","E_{T}(GeV)");
  fHistNPPElectronAccEtaE = CreateEtaEHisto2D("fHistNPPElectronAccEtaE_","MC non-primary electrons","#");
  fHistNPPElectronAccEtaPt = CreateEtaPtHisto2D("fHistNPPElectronAccEtaPt_","MC non-primary electrons","#");
	
  fHistNPPGammaAccEtaEET = CreateEtaEHisto2D("fHistNPPGammaAccEtaEET_","MC E_{T}, non-primary gammas","E_{T}(GeV)"); 
  fHistNPPGammaAccEtaPtET = CreateEtaPtHisto2D("fHistNPPGammaAccEtaPtET_","MC E_{T}, non-primary gammas","E_{T}(GeV)"); 
  fHistNPPGammaAccEtaE = CreateEtaEHisto2D("fHistNPPGammaAccEtaE_","MC non-primary gammas","#"); 
  fHistNPPGammaAccEtaPt = CreateEtaPtHisto2D("fHistNPPGammaAccEtaPt_","MC non-primary gammas","#"); 
	
  fHistElectronRecEtaEET = CreateEtaEHisto2D("fHistElectronRecEtaEET_","MC E_{T}, primary electrons","E_{T}(GeV)");
  fHistElectronRecEtaPtET = CreateEtaPtHisto2D("fHistElectronRecEtaPtET_","MC E_{T}, primary electrons","E_{T}(GeV)");
  fHistElectronRecEtaET = CreateEtaEtHisto2D("fHistElectronRecEtaET_","MC primary electrons","#");
  fHistElectronRecEtaE = CreateEtaEHisto2D("fHistElectronRecEtaE_","MC primary electrons","#");
  fHistElectronRecEtaPt = CreateEtaPtHisto2D("fHistElectronRecEtaPt_","MC primary electrons","#");
  histname = "fHistElectronRectotET_" + fHistogramNameSuffix;
  fHistElectronRectotET = new TH1F(histname.Data(),"total ET, MC primary electrons",fgNumOfEBins, fgEAxis);

  fHistConvElectronRecEtaEET = CreateEtaEHisto2D("fHistConvElectronRecEtaEET_","MC E_{T}, electrons from conversion","E_{T}(GeV)");
  fHistConvElectronRecEtaPtET = CreateEtaPtHisto2D("fHistConvElectronRecEtaPtET_","MC E_{T}, electrons from conversion","E_{T}(GeV)");
  fHistConvElectronRecEtaET = CreateEtaEtHisto2D("fHistConvElectronRecEtaET_","MC electrons from conversion","#");
  fHistConvElectronRecEtaE = CreateEtaEHisto2D("fHistConvElectronRecEtaE_","MC electrons from conversion","#");
  fHistConvElectronRecEtaPt = CreateEtaPtHisto2D("fHistConvElectronRecEtaPt_","MC electrons from conversion","#");
  histname = "fHistConvElectronRectotET_" + fHistogramNameSuffix;
  fHistConvElectronRectotET = new TH1F(histname.Data(),"total ET, MC electrons from conversion",fgNumOfEBins, fgEAxis);
	
  fHistScatElectronRecEtaEET = CreateEtaEHisto2D("fHistScatElectronRecEtaEET_","MC E_{T}, electrons from Scattering","E_{T}(GeV)");
  fHistScatElectronRecEtaPtET = CreateEtaPtHisto2D("fHistScatElectronRecEtaPtET_","MC E_{T}, electrons from Scattering","E_{T}(GeV)");
  fHistScatElectronRecEtaET = CreateEtaEtHisto2D("fHistScatElectronRecEtaET_","MC electrons from Scattering","#");
  fHistScatElectronRecEtaE = CreateEtaEHisto2D("fHistScatElectronRecEtaE_","MC electrons from Scattering","#");
  fHistScatElectronRecEtaPt = CreateEtaPtHisto2D("fHistScatElectronRecEtaPt_","MC electrons from Scattering","#");
  histname = "fHistScatElectronRectotET_" + fHistogramNameSuffix;
  fHistScatElectronRectotET = new TH1F(histname.Data(),"total ET, MC electrons from Scattering",fgNumOfEBins, fgEAxis);
	
  histname = "fHistTotElectronRectotET_" + fHistogramNameSuffix;
  fHistTotElectronRectotET = new TH1F(histname.Data(),"total ET, MC primary electrons",fgNumOfEBins, fgEAxis);

  fHistGammaRecEtaEET = CreateEtaEHisto2D("fHistGammaRecEtaEET_","MC E_{T}, primary gammas","E_{T}(GeV)"); 
  fHistGammaRecEtaPtET = CreateEtaPtHisto2D("fHistGammaRecEtaPtET_","MC E_{T}, primary gammas","E_{T}(GeV)"); 
  fHistGammaRecEtaET = CreateEtaEtHisto2D("fHistGammaRecEtaET_","MC primary gammas","#"); 
  fHistGammaRecEtaE = CreateEtaEHisto2D("fHistGammaRecEtaE_","MC primary gammas","#"); 
  fHistGammaRecEtaPt = CreateEtaPtHisto2D("fHistGammaRecEtaPt_","MC primary gammas","#"); 
  histname = "fHistGammaRectotET_" + fHistogramNameSuffix;
  fHistGammaRectotET = new TH1F(histname.Data(),"total ET, MC primary gammas",fgNumOfEBins, fgEAxis);

  fHistAnnihGammaRecEtaEET = CreateEtaEHisto2D("fHistAnnihGammaRecEtaEET_","MC E_{T}, Annihilation gammas","E_{T}(GeV)"); 
  fHistAnnihGammaRecEtaPtET = CreateEtaPtHisto2D("fHistAnnihGammaRecEtaPtET_","MC E_{T}, Annihilation gammas","E_{T}(GeV)"); 
  fHistAnnihGammaRecEtaET = CreateEtaEtHisto2D("fHistAnnihGammaRecEtaET_","MC Annihilation gammas","#"); 	
  fHistAnnihGammaRecEtaE = CreateEtaEHisto2D("fHistAnnihGammaRecEtaE_","MC Annihilation gammas","#"); 	
  fHistAnnihGammaRecEtaPt = CreateEtaPtHisto2D("fHistAnnihGammaRecEtaPt_","MC Annihilation gammas","#"); 
  histname = "fHistAnnihGammaRectotET_" + fHistogramNameSuffix;
  fHistAnnihGammaRectotET = new TH1F(histname.Data(),"total ET, MC Annihilation gammas",fgNumOfEBins, fgEAxis);
	
  fHistScatGammaRecEtaEET = CreateEtaEHisto2D("fHistScatGammaRecEtaEET_","MC E_{T}, Scattering gammas","E_{T}(GeV)"); 
  fHistScatGammaRecEtaPtET = CreateEtaPtHisto2D("fHistScatGammaRecEtaPtET_","MC E_{T}, Scattering gammas","E_{T}(GeV)"); 
  fHistScatGammaRecEtaET = CreateEtaEtHisto2D("fHistScatGammaRecEtaET_","MC Scattering gammas","#"); 	
  fHistScatGammaRecEtaE = CreateEtaEHisto2D("fHistScatGammaRecEtaE_","MC Scattering gammas","#"); 	
  fHistScatGammaRecEtaPt = CreateEtaPtHisto2D("fHistScatGammaRecEtaPt_","MC Scattering gammas","#"); 
  histname = "fHistScatGammaRectotET_" + fHistogramNameSuffix;
  fHistScatGammaRectotET = new TH1F(histname.Data(),"total ET, MC Scattering gammas",fgNumOfEBins, fgEAxis);
	
  histname = "fHistTotGammaRectotET_" + fHistogramNameSuffix;
  fHistTotGammaRectotET = new TH1F(histname.Data(),"total ET, MC primary gammas",fgNumOfEBins, fgEAxis);
	
  histname = "fHistTotEMRectotET_" + fHistogramNameSuffix;
  fHistTotEMRectotET = new TH1F(histname.Data(),"total electromagnetic ET",fgNumOfEBins, fgEAxis);

  fHistNPPElectronRecEtaEET = CreateEtaEHisto2D("fHistNPPElectronRecEtaEET_","MC E_{T}, non-primary electrons","E_{T}(GeV)");
  fHistNPPElectronRecEtaPtET = CreateEtaPtHisto2D("fHistNPPElectronRecEtaPtET_","MC E_{T}, non-primary electrons","E_{T}(GeV)");
  fHistNPPElectronRecEtaET = CreateEtaEtHisto2D("fHistNPPElectronRecEtaET_","MC non-primary electrons","#");
  fHistNPPElectronRecEtaE = CreateEtaEHisto2D("fHistNPPElectronRecEtaE_","MC non-primary electrons","#");
  fHistNPPElectronRecEtaPt = CreateEtaPtHisto2D("fHistNPPElectronRecEtaPt_","MC non-primary electrons","#");
  histname = "fHistNPPElectronRectotET_" + fHistogramNameSuffix;
  fHistNPPElectronRectotET = new TH1F(histname.Data(),"total ET, MC non-primary electrons",fgNumOfEBins, fgEAxis);
	
  fHistNPPGammaRecEtaEET = CreateEtaEHisto2D("fHistNPPGammaRecEtaEET_","MC E_{T}, non-primary gammas","E_{T}(GeV)"); 
  fHistNPPGammaRecEtaPtET = CreateEtaPtHisto2D("fHistNPPGammaRecEtaPtET_","MC E_{T}, non-primary gammas","E_{T}(GeV)"); 
  fHistNPPGammaRecEtaET = CreateEtaEtHisto2D("fHistNPPGammaRecEtaET_","MC non-primary gammas","#"); 
  fHistNPPGammaRecEtaE = CreateEtaEHisto2D("fHistNPPGammaRecEtaE_","MC non-primary gammas","#"); 
  fHistNPPGammaRecEtaPt = CreateEtaPtHisto2D("fHistNPPGammaRecEtaPt_","MC non-primary gammas","#"); 
  histname = "fHistNPPGammaRectotET_" + fHistogramNameSuffix;
  fHistNPPGammaRectotET = new TH1F(histname.Data(),"total ET, MC non-primary gammas",fgNumOfEBins, fgEAxis);

  histname = "fHistTotNPPEMRectotET_" + fHistogramNameSuffix;
  fHistTotNPPEMRectotET = new TH1F(histname.Data(),"total ET, MC non-primary electromagnetic",fgNumOfEBins, fgEAxis);

  fHistNPPPi0GammaRecEtaEET = CreateEtaEHisto2D("fHistNPPPi0GammaRecEtaEET_","MC E_{T}, non-primary gammas","E_{T}(GeV)"); 
  fHistNPPPi0GammaRecEtaPtET = CreateEtaPtHisto2D("fHistNPPPi0GammaRecEtaPtET_","MC E_{T}, non-primary gammas","E_{T}(GeV)"); 
  fHistNPPPi0GammaRecEtaET = CreateEtaEtHisto2D("fHistNPPPi0GammaRecEtaET_","MC non-primary gammas","#"); 
  fHistNPPPi0GammaRecEtaE = CreateEtaEHisto2D("fHistNPPPi0GammaRecEtaE_","MC non-primary gammas","#"); 
  fHistNPPPi0GammaRecEtaPt = CreateEtaPtHisto2D("fHistNPPPi0GammaRecEtaPt_","MC non-primary gammas","#"); 
  histname = "fHistNPPPi0GammaRectotET_" + fHistogramNameSuffix;
  fHistNPPPi0GammaRectotET = new TH1F(histname.Data(),"total ET, MC non-primary gammas",fgNumOfEBins, fgEAxis);
	
  fHistMuonEtaEET = CreateEtaEHisto2D("fHistMuonEtaEET_","MC E_{T}, primary Muons","E_{T}(GeV)");
  fHistMuonAccEtaEET = CreateEtaEHisto2D("fHistMuonAccEtaEET_","MC E_{T}, primary Muons, inside EMCal acceptance","E_{T}(GeV)");
  fHistMuonRecEtaEET = CreateEtaEHisto2D("fHistMuonRecEtaEET_","MC E_{T}, primary Muons, reconstructed","E_{T}(GeV)");
  fHistMuonMatchEtaEET = CreateEtaEHisto2D("fHistMuonMatchEtaEET_","MC E_{T}, primary Muons, tracking matched","E_{T}(GeV)");
	
  fHistMuonEtaPtET = CreateEtaPtHisto2D("fHistMuonEtaPtET_","MC E_{T}, primary Muons","E_{T}(GeV)");
  fHistMuonAccEtaPtET = CreateEtaPtHisto2D("fHistMuonAccEtaPtET_","MC E_{T}, primary Muons","E_{T}(GeV)");
  fHistMuonRecEtaPtET = CreateEtaPtHisto2D("fHistMuonRecEtaPtET_","MC E_{T}, primary Muons","E_{T}(GeV)");
  fHistMuonMatchEtaPtET = CreateEtaPtHisto2D("fHistMuonMatchEtaPtET_","MC E_{T}, primary Muons","E_{T}(GeV)");

  fHistMuonEtaET = CreateEtaEtHisto2D("fHistMuonEtaET_","MC primary Muons","#");
  fHistMuonAccEtaET = CreateEtaEtHisto2D("fHistMuonAccEtaET_","MC primary Muons","#");
  fHistMuonRecEtaET = CreateEtaEtHisto2D("fHistMuonRecEtaET_","MC primary Muons","#");
  fHistMuonMatchEtaET = CreateEtaEtHisto2D("fHistMuonMatchEtaET_","MC primary Muons","#");
	
  fHistMuonEtaE = CreateEtaEHisto2D("fHistMuonEtaE_","MC primary Muons","#");
  fHistMuonAccEtaE = CreateEtaEHisto2D("fHistMuonAccEtaE_","MC primary Muons","#");
  fHistMuonRecEtaE = CreateEtaEHisto2D("fHistMuonRecEtaE_","MC primary Muons","#");
  fHistMuonMatchEtaE = CreateEtaEHisto2D("fHistMuonMatchEtaE_","MC primary Muons","#");

  fHistMuonEtaPt = CreateEtaPtHisto2D("fHistMuonEtaPt_","MC primary Muons","#");
  fHistMuonAccEtaPt = CreateEtaPtHisto2D("fHistMuonAccEtaPt_","MC primary Muons","#");
  fHistMuonRecEtaPt = CreateEtaPtHisto2D("fHistMuonRecEtaPt_","MC primary Muons","#");
  fHistMuonMatchEtaPt = CreateEtaPtHisto2D("fHistMuonMatchEtaPt_","MC primary Muons","#");

  histname = "fHistMuontotET_" + fHistogramNameSuffix;
  fHistMuontotET = new TH1F(histname.Data(),"total ET, MC primary Muons",fgNumOfEBins, fgEAxis);
  histname = "fHistMuonAcctotET_" + fHistogramNameSuffix;
  fHistMuonAcctotET = new TH1F(histname.Data(),"total ET, MC primary Muons",fgNumOfEBins, fgEAxis);
  histname = "fHistMuonRectotET_" + fHistogramNameSuffix;
  fHistMuonRectotET = new TH1F(histname.Data(),"total ET, MC primary Muons",fgNumOfEBins, fgEAxis);
  histname = "fHistMuonMatchtotET_" + fHistogramNameSuffix;
  fHistMuonMatchtotET = new TH1F(histname.Data(),"total ET, MC primary Muons",fgNumOfEBins, fgEAxis);

  histname = "fHistMuonRectotETDep_" + fHistogramNameSuffix;
  fHistMuonRectotETDep = new TH1F(histname.Data(),"total ET, MC primary Muons",fgNumOfEBins, fgEAxis);
  histname = "fHistMuonMatchtotETDep_" + fHistogramNameSuffix;
  fHistMuonMatchtotETDep = new TH1F(histname.Data(),"total ET, MC primary Muons",fgNumOfEBins, fgEAxis);
		
  fHistMuonRecEtaEDepETDep = CreateEtaEHisto2D("fHistMuonRecEtaEDepETDep_","MC E_{T}, primary Muons, reconstructed","E_{T} dep (GeV)");
  fHistMuonMatchEtaEDepETDep = CreateEtaEHisto2D("fHistMuonMatchEtaEDepETDep_","MC E_{T}, primary Muons, tracking matched","E_{T} dep (GeV)");
	
  fHistMuonRecEtaPtETDep = CreateEtaPtHisto2D("fHistMuonRecEtaPtETDep_","MC E_{T}, primary Muons","E_{T} dep (GeV)");
  fHistMuonMatchEtaPtETDep = CreateEtaPtHisto2D("fHistMuonMatchEtaPtETDep_","MC E_{T}, primary Muons","E_{T} dep(GeV)");
	
  fHistMuonRecEtaETDep = CreateEtaEtHisto2D("fHistMuonRecEtaETDep_","MC primary Muons","#");
  fHistMuonMatchEtaETDep = CreateEtaEtHisto2D("fHistMuonMatchEtaETDep_","MC primary Muons","#");

  fHistMuonRecResEET = CreateResEHisto2D("fHistMuonRecResEET_","MC E_{T}, primary Muons","E_{T}(GeV)");
  fHistMuonRecResPtET = CreateResPtHisto2D("fHistMuonRecResPtET_","MC E_{T}, primary Muons","E_{T}(GeV)");
  fHistMuonRecResE = CreateResEHisto2D("fHistMuonRecResE_","MC primary Muons","#");
  fHistMuonRecResPt  = CreateResPtHisto2D("fHistMuonRecResPt_","MC primary Muons","#");
  fHistMuonRecResEDepETDep = CreateResEHisto2D("fHistMuonRecResEDepETDep_","MC E_{T}, primary Muons","E_{T} dep (GeV)");
  fHistMuonRecResPtETDep = CreateResPtHisto2D("fHistMuonRecResPtETDep_","MC E_{T}, primary Muons","E_{T} dep (GeV)");
	
  fHistPionEtaEET = CreateEtaEHisto2D("fHistPionEtaEET_","MC E_{T}, primary Pions","E_{T}(GeV)");
  fHistPionAccEtaEET = CreateEtaEHisto2D("fHistPionAccEtaEET_","MC E_{T}, primary Pions, inside EMCal acceptance","E_{T}(GeV)");
  fHistPionRecEtaEET = CreateEtaEHisto2D("fHistPionRecEtaEET_","MC E_{T}, primary Pions, reconstructed","E_{T}(GeV)");
  fHistPionMatchEtaEET = CreateEtaEHisto2D("fHistPionMatchEtaEET_","MC E_{T}, primary Pions, tracking matched","E_{T}(GeV)");
	
  fHistPionEtaPtET = CreateEtaPtHisto2D("fHistPionEtaPtET_","MC E_{T}, primary Pions","E_{T}(GeV)");
  fHistPionAccEtaPtET = CreateEtaPtHisto2D("fHistPionAccEtaPtET_","MC E_{T}, primary Pions","E_{T}(GeV)");
  fHistPionRecEtaPtET = CreateEtaPtHisto2D("fHistPionRecEtaPtET_","MC E_{T}, primary Pions","E_{T}(GeV)");
  fHistPionMatchEtaPtET = CreateEtaPtHisto2D("fHistPionMatchEtaPtET_","MC E_{T}, primary Pions","E_{T}(GeV)");
	
  fHistPionEtaET = CreateEtaEtHisto2D("fHistPionEtaET_","MC primary Pions","#");
  fHistPionAccEtaET = CreateEtaEtHisto2D("fHistPionAccEtaET_","MC primary Pions","#");
  fHistPionRecEtaET = CreateEtaEtHisto2D("fHistPionRecEtaET_","MC primary Pions","#");
  fHistPionMatchEtaET = CreateEtaEtHisto2D("fHistPionMatchEtaET_","MC primary Pions","#");
	
  fHistPionEtaE = CreateEtaEHisto2D("fHistPionEtaE_","MC primary Pions","#");
  fHistPionAccEtaE = CreateEtaEHisto2D("fHistPionAccEtaE_","MC primary Pions","#");
  fHistPionRecEtaE = CreateEtaEHisto2D("fHistPionRecEtaE_","MC primary Pions","#");
  fHistPionMatchEtaE = CreateEtaEHisto2D("fHistPionMatchEtaE_","MC primary Pions","#");
	
  fHistPionEtaPt = CreateEtaPtHisto2D("fHistPionEtaPt_","MC primary Pions","#");
  fHistPionAccEtaPt = CreateEtaPtHisto2D("fHistPionAccEtaPt_","MC primary Pions","#");
  fHistPionRecEtaPt = CreateEtaPtHisto2D("fHistPionRecEtaPt_","MC primary Pions","#");
  fHistPionMatchEtaPt = CreateEtaPtHisto2D("fHistPionMatchEtaPt_","MC primary Pions","#");

  histname = "fHistPiontotET_" + fHistogramNameSuffix;
  fHistPiontotET = new TH1F(histname.Data(),"total ET, MC primary Pions",fgNumOfEBins, fgEAxis);
  histname = "fHistPionAcctotET_" + fHistogramNameSuffix;
  fHistPionAcctotET = new TH1F(histname.Data(),"total ET, MC primary Pions",fgNumOfEBins, fgEAxis);
  histname = "fHistPionRectotET_" + fHistogramNameSuffix;
  fHistPionRectotET = new TH1F(histname.Data(),"total ET, MC primary Pions",fgNumOfEBins, fgEAxis);
  histname = "fHistPionMatchtotET_" + fHistogramNameSuffix;
  fHistPionMatchtotET = new TH1F(histname.Data(),"total ET, MC primary Pions",fgNumOfEBins, fgEAxis);
	
  histname = "fHistPionRectotETDep_" + fHistogramNameSuffix;
  fHistPionRectotETDep = new TH1F(histname.Data(),"total ET, MC primary Pions",fgNumOfEBins, fgEAxis);
  histname = "fHistPionMatchtotETDep_" + fHistogramNameSuffix;
  fHistPionMatchtotETDep = new TH1F(histname.Data(),"total ET, MC primary Pions",fgNumOfEBins, fgEAxis);
	
  fHistPionRecEtaEDepETDep = CreateEtaEHisto2D("fHistPionRecEtaEDepETDep_","MC E_{T}, primary Pions, reconstructed","E_{T} dep (GeV)");
  fHistPionMatchEtaEDepETDep = CreateEtaEHisto2D("fHistPionMatchEtaEDepETDep_","MC E_{T}, primary Pions, tracking matched","E_{T} dep (GeV)");
	
  fHistPionRecEtaPtETDep = CreateEtaPtHisto2D("fHistPionRecEtaPtETDep_","MC E_{T}, primary Pions","E_{T} dep (GeV)");
  fHistPionMatchEtaPtETDep = CreateEtaPtHisto2D("fHistPionMatchEtaPtETDep_","MC E_{T}, primary Pions","E_{T} dep(GeV)");
	
  fHistPionRecEtaETDep = CreateEtaEtHisto2D("fHistPionRecEtaETDep_","MC primary Pions","#");
  fHistPionMatchEtaETDep = CreateEtaEtHisto2D("fHistPionMatchEtaETDep_","MC primary Pions","#");
	
  fHistPionRecResEET = CreateResEHisto2D("fHistPionRecResEET_","MC E_{T}, primary Pions","E_{T}(GeV)");
  fHistPionRecResPtET = CreateResPtHisto2D("fHistPionRecResPtET_","MC E_{T}, primary Pions","E_{T}(GeV)");
  fHistPionRecResE = CreateResEHisto2D("fHistPionRecResE_","MC primary Pions","#");
  fHistPionRecResPt  = CreateResPtHisto2D("fHistPionRecResPt_","MC primary Pions","#");
  fHistPionRecResEDepETDep = CreateResEHisto2D("fHistPionRecResEDepETDep_","MC E_{T}, primary Pions","E_{T} dep (GeV)");
  fHistPionRecResPtETDep = CreateResPtHisto2D("fHistPionRecResPtETDep_","MC E_{T}, primary Pions","E_{T} dep (GeV)");
	
  fHistKaonEtaEET = CreateEtaEHisto2D("fHistKaonEtaEET_","MC E_{T}, primary Kaons","E_{T}(GeV)");
  fHistKaonAccEtaEET = CreateEtaEHisto2D("fHistKaonAccEtaEET_","MC E_{T}, primary Kaons, inside EMCal acceptance","E_{T}(GeV)");
  fHistKaonRecEtaEET = CreateEtaEHisto2D("fHistKaonRecEtaEET_","MC E_{T}, primary Kaons, reconstructed","E_{T}(GeV)");
  fHistKaonMatchEtaEET = CreateEtaEHisto2D("fHistKaonMatchEtaEET_","MC E_{T}, primary Kaons, tracking matched","E_{T}(GeV)");
	
  fHistKaonEtaPtET = CreateEtaPtHisto2D("fHistKaonEtaPtET_","MC E_{T}, primary Kaons","E_{T}(GeV)");
  fHistKaonAccEtaPtET = CreateEtaPtHisto2D("fHistKaonAccEtaPtET_","MC E_{T}, primary Kaons","E_{T}(GeV)");
  fHistKaonRecEtaPtET = CreateEtaPtHisto2D("fHistKaonRecEtaPtET_","MC E_{T}, primary Kaons","E_{T}(GeV)");
  fHistKaonMatchEtaPtET = CreateEtaPtHisto2D("fHistKaonMatchEtaPtET_","MC E_{T}, primary Kaons","E_{T}(GeV)");
	
  fHistKaonEtaET = CreateEtaEtHisto2D("fHistKaonEtaET_","MC primary Kaons","#");
  fHistKaonAccEtaET = CreateEtaEtHisto2D("fHistKaonAccEtaET_","MC primary Kaons","#");
  fHistKaonRecEtaET = CreateEtaEtHisto2D("fHistKaonRecEtaET_","MC primary Kaons","#");
  fHistKaonMatchEtaET = CreateEtaEtHisto2D("fHistKaonMatchEtaET_","MC primary Kaons","#");
	
  fHistKaonEtaE = CreateEtaEHisto2D("fHistKaonEtaE_","MC primary Kaons","#");
  fHistKaonAccEtaE = CreateEtaEHisto2D("fHistKaonAccEtaE_","MC primary Kaons","#");
  fHistKaonRecEtaE = CreateEtaEHisto2D("fHistKaonRecEtaE_","MC primary Kaons","#");
  fHistKaonMatchEtaE = CreateEtaEHisto2D("fHistKaonMatchEtaE_","MC primary Kaons","#");
	
  fHistKaonEtaPt = CreateEtaPtHisto2D("fHistKaonEtaPt_","MC primary Kaons","#");
  fHistKaonAccEtaPt = CreateEtaPtHisto2D("fHistKaonAccEtaPt_","MC primary Kaons","#");
  fHistKaonRecEtaPt = CreateEtaPtHisto2D("fHistKaonRecEtaPt_","MC primary Kaons","#");
  fHistKaonMatchEtaPt = CreateEtaPtHisto2D("fHistKaonMatchEtaPt_","MC primary Kaons","#");
	
  histname = "fHistKaontotET_" + fHistogramNameSuffix;
  fHistKaontotET = new TH1F(histname.Data(),"total ET, MC primary Kaons",fgNumOfEBins, fgEAxis);
  histname = "fHistKaonAcctotET_" + fHistogramNameSuffix;
  fHistKaonAcctotET = new TH1F(histname.Data(),"total ET, MC primary Kaons",fgNumOfEBins, fgEAxis);
  histname = "fHistKaonRectotET_" + fHistogramNameSuffix;
  fHistKaonRectotET = new TH1F(histname.Data(),"total ET, MC primary Kaons",fgNumOfEBins, fgEAxis);
  histname = "fHistKaonMatchtotET_" + fHistogramNameSuffix;
  fHistKaonMatchtotET = new TH1F(histname.Data(),"total ET, MC primary Kaons",fgNumOfEBins, fgEAxis);
	
  histname = "fHistKaonRectotETDep_" + fHistogramNameSuffix;
  fHistKaonRectotETDep = new TH1F(histname.Data(),"total ET, MC primary Kaons",fgNumOfEBins, fgEAxis);
  histname = "fHistKaonMatchtotETDep_" + fHistogramNameSuffix;
  fHistKaonMatchtotETDep = new TH1F(histname.Data(),"total ET, MC primary Kaons",fgNumOfEBins, fgEAxis);
	
  fHistKaonRecEtaEDepETDep = CreateEtaEHisto2D("fHistKaonRecEtaEDepETDep_","MC E_{T}, primary Kaons, reconstructed","E_{T} dep (GeV)");
  fHistKaonMatchEtaEDepETDep = CreateEtaEHisto2D("fHistKaonMatchEtaEDepETDep_","MC E_{T}, primary Kaons, tracking matched","E_{T} dep (GeV)");
	
  fHistKaonRecEtaPtETDep = CreateEtaPtHisto2D("fHistKaonRecEtaPtETDep_","MC E_{T}, primary Kaons","E_{T} dep (GeV)");
  fHistKaonMatchEtaPtETDep = CreateEtaPtHisto2D("fHistKaonMatchEtaPtETDep_","MC E_{T}, primary Kaons","E_{T} dep(GeV)");
	
  fHistKaonRecEtaETDep = CreateEtaEtHisto2D("fHistKaonRecEtaETDep_","MC primary Kaons","#");
  fHistKaonMatchEtaETDep = CreateEtaEtHisto2D("fHistKaonMatchEtaETDep_","MC primary Kaons","#");
	
  fHistKaonRecResEET = CreateResEHisto2D("fHistKaonRecResEET_","MC E_{T}, primary Kaons","E_{T}(GeV)");
  fHistKaonRecResPtET = CreateResPtHisto2D("fHistKaonRecResPtET_","MC E_{T}, primary Kaons","E_{T}(GeV)");
  fHistKaonRecResE = CreateResEHisto2D("fHistKaonRecResE_","MC primary Kaons","#");
  fHistKaonRecResPt  = CreateResPtHisto2D("fHistKaonRecResPt_","MC primary Kaons","#");	
  fHistKaonRecResEDepETDep = CreateResEHisto2D("fHistKaonRecResEDepETDep_","MC E_{T}, primary Kaons","E_{T} dep (GeV)");
  fHistKaonRecResPtETDep = CreateResPtHisto2D("fHistKaonRecResPtETDep_","MC E_{T}, primary Kaons","E_{T} dep (GeV)");

  fHistProtonEtaEET = CreateEtaEHisto2D("fHistProtonEtaEET_","MC E_{T}, primary Protons","E_{T}(GeV)");
  fHistProtonAccEtaEET = CreateEtaEHisto2D("fHistProtonAccEtaEET_","MC E_{T}, primary Protons, inside EMCal acceptance","E_{T}(GeV)");
  fHistProtonRecEtaEET = CreateEtaEHisto2D("fHistProtonRecEtaEET_","MC E_{T}, primary Protons, reconstructed","E_{T}(GeV)");
  fHistProtonMatchEtaEET = CreateEtaEHisto2D("fHistProtonMatchEtaEET_","MC E_{T}, primary Protons, tracking matched","E_{T}(GeV)");
	
  fHistProtonEtaPtET = CreateEtaPtHisto2D("fHistProtonEtaPtET_","MC E_{T}, primary Protons","E_{T}(GeV)");
  fHistProtonAccEtaPtET = CreateEtaPtHisto2D("fHistProtonAccEtaPtET_","MC E_{T}, primary Protons","E_{T}(GeV)");
  fHistProtonRecEtaPtET = CreateEtaPtHisto2D("fHistProtonRecEtaPtET_","MC E_{T}, primary Protons","E_{T}(GeV)");
  fHistProtonMatchEtaPtET = CreateEtaPtHisto2D("fHistProtonMatchEtaPtET_","MC E_{T}, primary Protons","E_{T}(GeV)");
	
  fHistProtonEtaET = CreateEtaEtHisto2D("fHistProtonEtaET_","MC primary Protons","#");
  fHistProtonAccEtaET = CreateEtaEtHisto2D("fHistProtonAccEtaET_","MC primary Protons","#");
  fHistProtonRecEtaET = CreateEtaEtHisto2D("fHistProtonRecEtaET_","MC primary Protons","#");
  fHistProtonMatchEtaET = CreateEtaEtHisto2D("fHistProtonMatchEtaET_","MC primary Protons","#");
	
  fHistProtonEtaE = CreateEtaEHisto2D("fHistProtonEtaE_","MC primary Protons","#");
  fHistProtonAccEtaE = CreateEtaEHisto2D("fHistProtonAccEtaE_","MC primary Protons","#");
  fHistProtonRecEtaE = CreateEtaEHisto2D("fHistProtonRecEtaE_","MC primary Protons","#");
  fHistProtonMatchEtaE = CreateEtaEHisto2D("fHistProtonMatchEtaE_","MC primary Protons","#");
	
  fHistProtonEtaPt = CreateEtaPtHisto2D("fHistProtonEtaPt_","MC primary Protons","#");
  fHistProtonAccEtaPt = CreateEtaPtHisto2D("fHistProtonAccEtaPt_","MC primary Protons","#");
  fHistProtonRecEtaPt = CreateEtaPtHisto2D("fHistProtonRecEtaPt_","MC primary Protons","#");
  fHistProtonMatchEtaPt = CreateEtaPtHisto2D("fHistProtonMatchEtaPt_","MC primary Protons","#");
	
  histname = "fHistProtontotET_" + fHistogramNameSuffix;
  fHistProtontotET = new TH1F(histname.Data(),"total ET, MC primary Protons",fgNumOfEBins, fgEAxis);
  histname = "fHistProtonAcctotET_" + fHistogramNameSuffix;
  fHistProtonAcctotET = new TH1F(histname.Data(),"total ET, MC primary Protons",fgNumOfEBins, fgEAxis);
  histname = "fHistProtonRectotET_" + fHistogramNameSuffix;
  fHistProtonRectotET = new TH1F(histname.Data(),"total ET, MC primary Protons",fgNumOfEBins, fgEAxis);
  histname = "fHistProtonMatchtotET_" + fHistogramNameSuffix;
  fHistProtonMatchtotET = new TH1F(histname.Data(),"total ET, MC primary Protons",fgNumOfEBins, fgEAxis);
	
  histname = "fHistProtonRectotETDep_" + fHistogramNameSuffix;
  fHistProtonRectotETDep = new TH1F(histname.Data(),"total ET, MC primary Protons",fgNumOfEBins, fgEAxis);
  histname = "fHistProtonMatchtotETDep_" + fHistogramNameSuffix;
  fHistProtonMatchtotETDep = new TH1F(histname.Data(),"total ET, MC primary Protons",fgNumOfEBins, fgEAxis);
	
  fHistProtonRecEtaEDepETDep = CreateEtaEHisto2D("fHistProtonRecEtaEDepETDep_","MC E_{T}, primary Protons, reconstructed","E_{T} dep (GeV)");
  fHistProtonMatchEtaEDepETDep = CreateEtaEHisto2D("fHistProtonMatchEtaEDepETDep_","MC E_{T}, primary Protons, tracking matched","E_{T} dep (GeV)");
	
  fHistProtonRecEtaPtETDep = CreateEtaPtHisto2D("fHistProtonRecEtaPtETDep_","MC E_{T}, primary Protons","E_{T} dep (GeV)");
  fHistProtonMatchEtaPtETDep = CreateEtaPtHisto2D("fHistProtonMatchEtaPtETDep_","MC E_{T}, primary Protons","E_{T} dep(GeV)");
	
  fHistProtonRecEtaETDep = CreateEtaEtHisto2D("fHistProtonRecEtaETDep_","MC primary Protons","#");
  fHistProtonMatchEtaETDep = CreateEtaEtHisto2D("fHistProtonMatchEtaETDep_","MC primary Protons","#");
	
  fHistProtonRecResEET = CreateResEHisto2D("fHistProtonRecResEET_","MC E_{T}, primary Protons","E_{T}(GeV)");
  fHistProtonRecResPtET = CreateResPtHisto2D("fHistProtonRecResPtET_","MC E_{T}, primary Protons","E_{T}(GeV)");
  fHistProtonRecResE = CreateResEHisto2D("fHistProtonRecResE_","MC primary Protons","#");
  fHistProtonRecResPt  = CreateResPtHisto2D("fHistProtonRecResPt_","MC primary Protons","#");
  fHistProtonRecResEDepETDep = CreateResEHisto2D("fHistProtonRecResEDepETDep_","MC E_{T}, primary Protons","E_{T} dep (GeV)");
  fHistProtonRecResPtETDep = CreateResPtHisto2D("fHistProtonRecResPtETDep_","MC E_{T}, primary Protons","E_{T} dep (GeV)");
	
  histname = "fHistTotChargedtotET_" + fHistogramNameSuffix;
  fHistTotChargedtotET = new TH1F(histname.Data(),"total ET, MC primary charged particles",fgNumOfEBins, fgEAxis);
  histname = "fHistTotChargedAcctotET_" + fHistogramNameSuffix;
  fHistTotChargedAcctotET = new TH1F(histname.Data(),"total ET, MC primary charged particles",fgNumOfEBins, fgEAxis);
  histname = "fHistTotChargedRectotET_" + fHistogramNameSuffix;
  fHistTotChargedRectotET = new TH1F(histname.Data(),"total ET, MC primary charged particles",fgNumOfEBins, fgEAxis);
  histname = "fHistTotChargedMatchtotET_" + fHistogramNameSuffix;
  fHistTotChargedMatchtotET = new TH1F(histname.Data(),"total ET, MC primary charged particles",fgNumOfEBins, fgEAxis);
	
  histname = "fHistTotChargedRectotETDep_" + fHistogramNameSuffix;
  fHistTotChargedRectotETDep = new TH1F(histname.Data(),"total ET, MC primary charged particles",fgNumOfEBins, fgEAxis);
  histname = "fHistTotChargedMatchtotETDep_" + fHistogramNameSuffix;
  fHistTotChargedMatchtotETDep = new TH1F(histname.Data(),"total ET, MC primary charged particles",fgNumOfEBins, fgEAxis);
	
  fHistNeutronEtaEET = CreateEtaEHisto2D("fHistNeutronEtaEET_","MC E_{T}, primary Neutrons","E_{T}(GeV)");
  fHistNeutronAccEtaEET = CreateEtaEHisto2D("fHistNeutronAccEtaEET_","MC E_{T}, primary Neutrons, inside EMCal acceptance","E_{T}(GeV)");
  fHistNeutronRecEtaEET = CreateEtaEHisto2D("fHistNeutronRecEtaEET_","MC E_{T}, primary Neutrons, reconstructed","E_{T}(GeV)");
	
  fHistNeutronEtaPtET = CreateEtaPtHisto2D("fHistNeutronEtaPtET_","MC E_{T}, primary Neutrons","E_{T}(GeV)");
  fHistNeutronAccEtaPtET = CreateEtaPtHisto2D("fHistNeutronAccEtaPtET_","MC E_{T}, primary Neutrons","E_{T}(GeV)");
  fHistNeutronRecEtaPtET = CreateEtaPtHisto2D("fHistNeutronRecEtaPtET_","MC E_{T}, primary Neutrons","E_{T}(GeV)");
	
  fHistNeutronEtaET = CreateEtaEtHisto2D("fHistNeutronEtaET_","MC primary Neutrons","#");
  fHistNeutronAccEtaET = CreateEtaEtHisto2D("fHistNeutronAccEtaET_","MC primary Neutrons","#");
  fHistNeutronRecEtaET = CreateEtaEtHisto2D("fHistNeutronRecEtaET_","MC primary Neutrons","#");
	
  fHistNeutronEtaE = CreateEtaEHisto2D("fHistNeutronEtaE_","MC primary Neutrons","#");
  fHistNeutronAccEtaE = CreateEtaEHisto2D("fHistNeutronAccEtaE_","MC primary Neutrons","#");
  fHistNeutronRecEtaE = CreateEtaEHisto2D("fHistNeutronRecEtaE_","MC primary Neutrons","#");
	
  fHistNeutronEtaPt = CreateEtaPtHisto2D("fHistNeutronEtaPt_","MC primary Neutrons","#");
  fHistNeutronAccEtaPt = CreateEtaPtHisto2D("fHistNeutronAccEtaPt_","MC primary Neutrons","#");
  fHistNeutronRecEtaPt = CreateEtaPtHisto2D("fHistNeutronRecEtaPt_","MC primary Neutrons","#");
	
  histname = "fHistNeutrontotET_" + fHistogramNameSuffix;
  fHistNeutrontotET = new TH1F(histname.Data(),"total ET, MC primary Neutrons",fgNumOfEBins, fgEAxis);
  histname = "fHistNeutronAcctotET_" + fHistogramNameSuffix;
  fHistNeutronAcctotET = new TH1F(histname.Data(),"total ET, MC primary Neutrons",fgNumOfEBins, fgEAxis);
  histname = "fHistNeutronRectotET_" + fHistogramNameSuffix;
  fHistNeutronRectotET = new TH1F(histname.Data(),"total ET, MC primary Neutrons",fgNumOfEBins, fgEAxis);
  histname = "fHistNeutronRectotETDep_" + fHistogramNameSuffix;
  fHistNeutronRectotETDep = new TH1F(histname.Data(),"total ET, MC primary Neutrons",fgNumOfEBins, fgEAxis);	
	
  fHistNeutronRecEtaEDepETDep = CreateEtaEHisto2D("fHistNeutronRecEtaEDepETDep_","MC E_{T}, primary Neutrons, reconstructed","E_{T} dep (GeV)");
  fHistNeutronRecEtaETDep = CreateEtaEtHisto2D("fHistNeutronRecEtaETDep_","MC primary Neutrons","#");
	
  fHistNeutronRecEtaPtETDep = CreateEtaPtHisto2D("fHistNeutronRecEtaPtETDep_","MC E_{T}, primary Neutrons","E_{T} dep (GeV)");
	
  fHistK0EtaEET = CreateEtaEHisto2D("fHistK0EtaEET_","MC E_{T}, K0S daughters","E_{T}(GeV)");
  fHistK0RecEtaEET = CreateEtaEHisto2D("fHistK0RecEtaEET_","MC E_{T}, K0S daughters, reconstructed","E_{T}(GeV)");
	
  fHistK0EtaPtET = CreateEtaPtHisto2D("fHistK0EtaPtET_","MC E_{T}, K0S daughters","E_{T}(GeV)");
  fHistK0RecEtaPtET = CreateEtaPtHisto2D("fHistK0RecEtaPtET_","MC E_{T}, K0S daughters","E_{T}(GeV)");
	
  fHistK0EtaET = CreateEtaEtHisto2D("fHistK0EtaET_","MC K0S daughters","#");
  fHistK0RecEtaET = CreateEtaEtHisto2D("fHistK0RecEtaET_","MC K0S daughters","#");
	
  fHistK0EtaE = CreateEtaEHisto2D("fHistK0EtaE_","MC K0S daughters","#");
  fHistK0RecEtaE = CreateEtaEHisto2D("fHistK0RecEtaE_","MC K0S daughters","#");
	
  fHistK0EtaPt = CreateEtaPtHisto2D("fHistK0EtaPt_","MC K0S daughters","#");
  fHistK0RecEtaPt = CreateEtaPtHisto2D("fHistK0RecEtaPt_","MC K0S daughters","#");
	
  histname = "fHistK0totET_" + fHistogramNameSuffix;
  fHistK0totET = new TH1F(histname.Data(),"total ET, MC K0s daughters",fgNumOfEBins, fgEAxis);
  histname = "fHistK0RectotET_" + fHistogramNameSuffix;
  fHistK0RectotET = new TH1F(histname.Data(),"total ET, MC K0s daughters",fgNumOfEBins, fgEAxis);	
  histname = "fHistK0RectotETDep_" + fHistogramNameSuffix;
  fHistK0RectotETDep = new TH1F(histname.Data(),"total ET, MC K0s daughters",fgNumOfEBins, fgEAxis);	
	
  fHistK0RecEtaEDepETDep = CreateEtaEHisto2D("fHistK0RecEtaEDepETDep_","MC E_{T}, MC K0s daughters, reconstructed","E_{T} dep (GeV)");
  fHistK0RecEtaETDep = CreateEtaEtHisto2D("fHistK0RecEtaETDep_","MC K0s daughters","#");
	
  fHistK0RecEtaPtETDep = CreateEtaPtHisto2D("fHistK0RecEtaPtETDep_","MC E_{T}, MC K0s daughters","E_{T} dep (GeV)");
	
  fHistLambdaEtaEET = CreateEtaEHisto2D("fHistLambdaEtaEET_","MC E_{T}, Lambda daughters","E_{T}(GeV)");
  fHistLambdaRecEtaEET = CreateEtaEHisto2D("fHistLambdaRecEtaEET_","MC E_{T}, Lambda daughters, reconstructed","E_{T}(GeV)");
	
  fHistLambdaEtaPtET = CreateEtaPtHisto2D("fHistLambdaEtaPtET_","MC E_{T}, Lambda daughters","E_{T}(GeV)");
  fHistLambdaRecEtaPtET = CreateEtaPtHisto2D("fHistLambdaRecEtaPtET_","MC E_{T}, Lambda daughters","E_{T}(GeV)");
	
  fHistLambdaEtaET = CreateEtaEtHisto2D("fHistLambdaEtaET_","MC Lambda daughters","#");
  fHistLambdaRecEtaET = CreateEtaEtHisto2D("fHistLambdaRecEtaET_","MC Lambda daughters","#");
	
  fHistLambdaEtaE = CreateEtaEHisto2D("fHistLambdaEtaE_","MC Lambda daughters","#");
  fHistLambdaRecEtaE = CreateEtaEHisto2D("fHistLambdaRecEtaE_","MC Lambda daughters","#");
	
  fHistLambdaEtaPt = CreateEtaPtHisto2D("fHistLambdaEtaPt_","MC Lambda daughters","#");
  fHistLambdaRecEtaPt = CreateEtaPtHisto2D("fHistLambdaRecEtaPt_","MC Lambda daughters","#");
		
  histname = "fHistLambdatotET_" + fHistogramNameSuffix;
  fHistLambdatotET = new TH1F(histname.Data(),"total ET, MC Lambdas daughters",fgNumOfEBins, fgEAxis);
  histname = "fHistLambdaRectotET_" + fHistogramNameSuffix;
  fHistLambdaRectotET = new TH1F(histname.Data(),"total ET, MC Lambdas daughters",fgNumOfEBins, fgEAxis);
  histname = "fHistLambdaRectotETDep_" + fHistogramNameSuffix;
  fHistLambdaRectotETDep = new TH1F(histname.Data(),"total ET, MC Lambdas daughters",fgNumOfEBins, fgEAxis);	
	
  fHistLambdaRecEtaEDepETDep = CreateEtaEHisto2D("fHistLambdaRecEtaEDepETDep_","MC E_{T}, MC Lambdas daughters, reconstructed","E_{T} dep (GeV)");
  fHistLambdaRecEtaETDep = CreateEtaEtHisto2D("fHistLambdaRecEtaETDep_","MC Lambdas daughters","#");
	
  fHistLambdaRecEtaPtETDep = CreateEtaPtHisto2D("fHistLambdaRecEtaPtETDep_","MC E_{T}, MC Lambdas daughters","E_{T} dep (GeV)");

  histname = "fHistTotNeutraltotET_" + fHistogramNameSuffix;
  fHistTotNeutraltotET = new TH1F(histname.Data(),"total ET, MC Lambdas daughters",fgNumOfEBins, fgEAxis);
  histname = "fHistTotNeutralRectotET_" + fHistogramNameSuffix;
  fHistTotNeutralRectotET = new TH1F(histname.Data(),"total ET, MC Lambdas daughters",fgNumOfEBins, fgEAxis);
  histname = "fHistTotNeutralRectotETDep_" + fHistogramNameSuffix;
  fHistTotNeutralRectotETDep = new TH1F(histname.Data(),"total ET, MC Lambdas daughters",fgNumOfEBins, fgEAxis);	
	
  histname = "fHistTotaltotET_" + fHistogramNameSuffix;
  fHistTotaltotET = new TH1F(histname.Data(),"total ET, all particles",fgNumOfEBins, fgEAxis);
  histname = "fHistTotalAcctotET_" + fHistogramNameSuffix;
  fHistTotalAcctotET = new TH1F(histname.Data(),"total ET, all particles",fgNumOfEBins, fgEAxis);
  histname = "fHistTotalRectotET_" + fHistogramNameSuffix;
  fHistTotalRectotET = new TH1F(histname.Data(),"total ET, all particles",fgNumOfEBins, fgEAxis);
  histname = "fHistTotalRectotETDep_" + fHistogramNameSuffix;
  fHistTotalRectotETDep = new TH1F(histname.Data(),"total ET, all particles",fgNumOfEBins, fgEAxis);	
	
  histname = "fHistAll_ERecvsMC_" + fHistogramNameSuffix;
  fHistAllERecEMC = new TH2F(histname.Data(),"E cluster Rec vs MC, all particles",fgNumOfEBins, fgEAxis,fgNumOfEBins, fgEAxis);
  fHistAllERecEMC->SetXTitle("E_{MC}(GeV)");
  fHistAllERecEMC->SetYTitle("E_{Rec}(GeV)");
	
  histname = "fHistElectron_ERecvsMC_" + fHistogramNameSuffix;
  fHistElectronERecEMC = new TH2F(histname.Data(),"E cluster Rec vs MC, Electrons",fgNumOfEBins, fgEAxis,fgNumOfEBins, fgEAxis);
  fHistElectronERecEMC->SetXTitle("E_{MC}(GeV)");
  fHistElectronERecEMC->SetYTitle("E_{Rec}(GeV)");
	
  histname = "fHistGamma_ERecvsMC_" + fHistogramNameSuffix;
  fHistGammaERecEMC = new TH2F(histname.Data(),"E cluster Rec vs MC, Gammas",fgNumOfEBins, fgEAxis,fgNumOfEBins, fgEAxis);
  fHistGammaERecEMC->SetXTitle("E_{MC}(GeV)");
  fHistGammaERecEMC->SetYTitle("E_{Rec}(GeV)");
	
  histname = "fHistAllPtRecPtMC_" + fHistogramNameSuffix;
  fHistAllPtRecPtMC = new TH2F(histname.Data(),"pt track Rec vs MC, all particles",fgNumOfEBins, fgEAxis,fgNumOfEBins, fgEAxis);
  fHistAllPtRecPtMC->SetXTitle("p_{T}^{MC}(GeV/c)");
  fHistAllPtRecPtMC->SetYTitle("p_{T}^{Rec}(GeV/c)");	
	
  histname = "fHistChargedRes_" + fHistogramNameSuffix;
  fHistChargedRes = new TH2F(histname.Data(),"#Delta#phi vs #Delta#eta (track projection - cluster position), charged particles",200,-0.1,0.1,200,-0.1,0.1);
  fHistChargedRes->SetXTitle("#Delta#phi");
  fHistChargedRes->SetYTitle("#Delta#eta");
	
  histname = "fHistChargedRes2_" + fHistogramNameSuffix;
  fHistChargedRes2 = new TH2F(histname.Data(),"#Delta#phi vs #Delta#eta (track projection - cluster position), charged particles",200,-0.1,0.1,200,-0.1,0.1);
  fHistChargedRes2->SetXTitle("#Delta#phi");
  fHistChargedRes2->SetYTitle("#Delta#eta");
	
  histname = "fHistChargedRes3_" + fHistogramNameSuffix;
  fHistChargedRes3 = new TH2F(histname.Data(),"#Delta#phi vs #Delta#eta (track projection - cluster position), charged particles",200,-0.1,0.1,200,-0.1,0.1);
  fHistChargedRes3->SetXTitle("#Delta#phi");
  fHistChargedRes3->SetYTitle("#Delta#eta");
	
  histname = "fHistNeutralRes_" + fHistogramNameSuffix;
  fHistNeutralRes = new TH2F(histname.Data(),"#Delta#phi vs #Delta#eta (track projection - cluster position), neutral particles",200,-0.1,0.1,200,-0.1,0.1);
  fHistNeutralRes->SetXTitle("#Delta#phi");
  fHistNeutralRes->SetYTitle("#Delta#eta");
	
  histname = "fHistElectronRes_" + fHistogramNameSuffix;
  fHistElectronRes = new TH2F(histname.Data(),"#Delta#phi vs #Delta#eta (track projection - cluster position, Electrons",200,-0.1,0.1,200,-0.1,0.1);
  fHistElectronRes->SetXTitle("#Delta#phi");
  fHistElectronRes->SetYTitle("#Delta#eta");
	
  histname = "fHistGammaRes_" + fHistogramNameSuffix;
  fHistGammaRes = new TH2F(histname.Data(),"#Delta#phi vs #Delta#eta (track projection - cluster position, Gammas",200,-0.1,0.1,200,-0.1,0.1);
  fHistGammaRes->SetXTitle("#Delta#phi");
  fHistGammaRes->SetYTitle("#Delta#eta");
	
  histname = "fHistIsInAcc_" + fHistogramNameSuffix;
  //fHistIsInAcc = new TH2F(histname.Data(),"X,Y position of particle projection inside EMCal",1201,-600.5,600.5,1201,-600.5,600.5);
  //fHistIsInAcc->SetXTitle("X (cm)");
  //fHistIsInAcc->SetYTitle("Y (cm)");
  fHistIsInAcc = new TH2F(histname.Data(),"#phhi, #eta position of particle projection inside EMCal",360,0.,360.,200,-1.,1.);
  fHistIsInAcc->SetXTitle("#phi");
  fHistIsInAcc->SetYTitle("#eta");
	
  histname = "fHistElectronFirstMother_" + fHistogramNameSuffix;
  fHistElectronFirstMother = new TH1F(histname.Data(),"Electron First Mother PDG Code",1201,-600.5,600.5);
  histname = "fHistElectronFirstMotherXY_" + fHistogramNameSuffix;
  fHistElectronFirstMotherXY = new TH2F(histname.Data(),"Electron Mother X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistElectronNDaughters_" + fHistogramNameSuffix;
  fHistElectronNDaughters = new TH1F(histname.Data(),"Number of Electron Daugthers, inside EMCal acceptance",11,-0.5,10.5);
  histname = "fHistElectronDaughters_" + fHistogramNameSuffix;
  fHistElectronDaughters = new TH1F(histname.Data(),"Electron Daugther PDG Code",1201,-600.5,600.5);
  histname = "fHistElectronDaughtersXY_" + fHistogramNameSuffix;
  fHistElectronDaughtersXY = new TH2F(histname.Data(),"Electron Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
	
  histname = "fHistElectronFirstMotherAcc_" + fHistogramNameSuffix;
  fHistElectronFirstMotherAcc = new TH1F(histname.Data(),"Electron First Mother PDG Code, inside EMCal acceptance",1201,-600.5,600.5);
  histname = "fHistElectronFirstMotherXYAcc_" + fHistogramNameSuffix;
  fHistElectronFirstMotherXYAcc = new TH2F(histname.Data(),"Electron Mother X,Y vertex position, inside EMCal acceptance",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistElectronNDaughtersAcc_" + fHistogramNameSuffix;
  fHistElectronNDaughtersAcc = new TH1F(histname.Data(),"Number of Electron Daugthers, inside EMCal acceptance",11,-0.5,10.5);
  histname = "fHistElectronDaughtersAcc_" + fHistogramNameSuffix;
  fHistElectronDaughtersAcc = new TH1F(histname.Data(),"Electron Daugther PDG Code, inside EMCal acceptance",1201,-600.5,600.5);
  histname = "fHistElectronDaughtersXYAcc_" + fHistogramNameSuffix;
  fHistElectronDaughtersXYAcc = new TH2F(histname.Data(),"Electron Daugther X,Y vertex position, inside EMCal acceptance",1201,-600.5,600.5,1201,-600.5,600.5);
	
  histname = "fHistElectronFirstMotherRec_" + fHistogramNameSuffix;
  fHistElectronFirstMotherRec = new TH1F(histname.Data(),"Reconstructed Electron First Mother PDG Code",1201,-600.5,600.5);
  histname = "fHistElectronFirstMotherXYRec_" + fHistogramNameSuffix;
  fHistElectronFirstMotherXYRec = new TH2F(histname.Data(),"Electron Mother X,Y vertex position, inside EMCal acceptance",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistElectronNDaughtersRec_" + fHistogramNameSuffix;
  fHistElectronNDaughtersRec = new TH1F(histname.Data(),"Number of Electron Daugthers, inside EMCal acceptance",11,-0.5,10.5);
  histname = "fHistElectronDaughtersRec_" + fHistogramNameSuffix;
  fHistElectronDaughtersRec = new TH1F(histname.Data(),"Electron Daugther PDG Code",1201,-600.5,600.5);
  histname = "fHistElectronDaughtersXYRec_" + fHistogramNameSuffix;
  fHistElectronDaughtersXYRec = new TH2F(histname.Data(),"Electron Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);

  histname = "fHistNPPElectronFirstMother_" + fHistogramNameSuffix;
  fHistNPPElectronFirstMother = new TH1F(histname.Data(),"Electron First Mother PDG Code",1201,-600.5,600.5);
  histname = "fHistNPPElectronFirstMotherXY_" + fHistogramNameSuffix;
  fHistNPPElectronFirstMotherXY = new TH2F(histname.Data(),"Electron Mother X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistNPPElectronNDaughters_" + fHistogramNameSuffix;
  fHistNPPElectronNDaughters = new TH1F(histname.Data(),"Number of Electron Daugthers",11,-0.5,10.5);
  histname = "fHistNPPElectronDaughters_" + fHistogramNameSuffix;
  fHistNPPElectronDaughters = new TH1F(histname.Data(),"Electron Daugther PDG Code",1201,-600.5,600.5);
  histname = "fHistNPPElectronDaughtersXY_" + fHistogramNameSuffix;
  fHistNPPElectronDaughtersXY = new TH2F(histname.Data(),"Electron Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
	
  histname = "fHistNPPElectronFirstMotherAcc_" + fHistogramNameSuffix;
  fHistNPPElectronFirstMotherAcc = new TH1F(histname.Data(),"Electron First Mother PDG Code, inside EMCal acceptance",1201,-600.5,600.5);
  histname = "fHistNPPElectronFirstMotherXYAcc_" + fHistogramNameSuffix;
  fHistNPPElectronFirstMotherXYAcc = new TH2F(histname.Data(),"Electron Mother X,Y vertex position, inside EMCal acceptance",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistNPPElectronNDaughtersAcc_" + fHistogramNameSuffix;
  fHistNPPElectronNDaughtersAcc = new TH1F(histname.Data(),"Number of Electron Daugthers, inside EMCal acceptance",11,-0.5,10.5);
  histname = "fHistNPPElectronDaughtersAcc_" + fHistogramNameSuffix;
  fHistNPPElectronDaughtersAcc = new TH1F(histname.Data(),"Electron Daugther PDG Code, inside EMCal acceptance",1201,-600.5,600.5);
  histname = "fHistNPPElectronDaughtersXYAcc_" + fHistogramNameSuffix;
  fHistNPPElectronDaughtersXYAcc = new TH2F(histname.Data(),"Electron Daugther X,Y vertex position, inside EMCal acceptance",1201,-600.5,600.5,1201,-600.5,600.5);
	
  histname = "fHistNPPElectronFirstMotherRec_" + fHistogramNameSuffix;
  fHistNPPElectronFirstMotherRec = new TH1F(histname.Data(),"Reconstructed Electron First Mother PDG Code",1201,-600.5,600.5);
  histname = "fHistNPPElectronFirstMotherXYRec_" + fHistogramNameSuffix;
  fHistNPPElectronFirstMotherXYRec = new TH2F(histname.Data(),"Electron Mother X,Y vertex position, inside EMCal acceptance",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistNPPElectronNDaughtersRec_" + fHistogramNameSuffix;
  fHistNPPElectronNDaughtersRec = new TH1F(histname.Data(),"Number of Electron Daugthers, inside EMCal acceptance",11,-0.5,10.5);
  histname = "fHistNPPElectronDaughtersRec_" + fHistogramNameSuffix;
  fHistNPPElectronDaughtersRec = new TH1F(histname.Data(),"Electron Daugther PDG Code",1201,-600.5,600.5);
  histname = "fHistNPPElectronDaughtersXYRec_" + fHistogramNameSuffix;
  fHistNPPElectronDaughtersXYRec = new TH2F(histname.Data(),"Electron Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
	
  histname = "fHistGammaFirstMother_" + fHistogramNameSuffix;
  fHistGammaFirstMother = new TH1F(histname.Data(),"Gamma First Mother PDG Code",1201,-600.5,600.5);
  histname = "fHistGammaFirstMotherXY_" + fHistogramNameSuffix;
  fHistGammaFirstMotherXY = new TH2F(histname.Data(),"Gamma Mother X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistGammaNDaughters_" + fHistogramNameSuffix;
  fHistGammaNDaughters = new TH1F(histname.Data(),"Number of Gamma Daugthers, inside EMCal acceptance",11,-0.5,10.5);
  histname = "fHistGammaDaughters_" + fHistogramNameSuffix;
  fHistGammaDaughters = new TH1F(histname.Data(),"Gamma Daugther PDG Code",1201,-600.5,600.5);
  histname = "fHistGammaDaughtersXY_" + fHistogramNameSuffix;
  fHistGammaDaughtersXY = new TH2F(histname.Data(),"Gamma Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistConvGammaDaughtersXY_" + fHistogramNameSuffix;
  fHistConvGammaDaughtersXY = new TH2F(histname.Data(),"Gamma Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistNonConvGammaDaughtersXY_" + fHistogramNameSuffix;
  fHistNonConvGammaDaughtersXY = new TH2F(histname.Data(),"Gamma Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
	
  histname = "fHistGammaFirstMotherAcc_" + fHistogramNameSuffix;
  fHistGammaFirstMotherAcc = new TH1F(histname.Data(),"Gamma First Mother PDG Code, inside EMCal acceptance",1201,-600.5,600.5);
  histname = "fHistGammaFirstMotherXYAcc_" + fHistogramNameSuffix;
  fHistGammaFirstMotherXYAcc = new TH2F(histname.Data(),"Gamma Mother X,Y vertex position, inside EMCal acceptance",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistGammaNDaughtersAcc_" + fHistogramNameSuffix;
  fHistGammaNDaughtersAcc = new TH1F(histname.Data(),"Number of Gamma Daugthers, inside EMCal acceptance",11,-0.5,10.5);
  histname = "fHistGammaDaughtersAcc_" + fHistogramNameSuffix;
  fHistGammaDaughtersAcc = new TH1F(histname.Data(),"Gamma Daugther PDG Code, inside EMCal acceptance",1201,-600.5,600.5);
  histname = "fHistGammaDaughtersXYAcc_" + fHistogramNameSuffix;
  fHistGammaDaughtersXYAcc = new TH2F(histname.Data(),"Gamma Daugther X,Y vertex position, inside EMCal acceptance",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistConvGammaDaughtersXYAcc_" + fHistogramNameSuffix;
  fHistConvGammaDaughtersXYAcc = new TH2F(histname.Data(),"Gamma Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistNonConvGammaDaughtersXYAcc_" + fHistogramNameSuffix;
  fHistNonConvGammaDaughtersXYAcc = new TH2F(histname.Data(),"Gamma Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
	
  histname = "fHistGammaFirstMotherRec_" + fHistogramNameSuffix;
  fHistGammaFirstMotherRec = new TH1F(histname.Data(),"Reconstructed Gamma First Mother PDG Code",1201,-600.5,600.5);
  histname = "fHistGammaFirstMotherXYRec_" + fHistogramNameSuffix;
  fHistGammaFirstMotherXYRec = new TH2F(histname.Data(),"Gamma Mother X,Y vertex position, inside EMCal acceptance",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistGammaNDaughtersRec_" + fHistogramNameSuffix;
  fHistGammaNDaughtersRec = new TH1F(histname.Data(),"Number of Gamma Daugthers, inside EMCal acceptance",11,-0.5,10.5);
  histname = "fHistGammaDaughtersRec_" + fHistogramNameSuffix;
  fHistGammaDaughtersRec = new TH1F(histname.Data(),"Gamma Daugther PDG Code",1201,-600.5,600.5);
  histname = "fHistGammaDaughtersXYRec_" + fHistogramNameSuffix;
  fHistGammaDaughtersXYRec = new TH2F(histname.Data(),"Gamma Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistConvGammaDaughtersXYRec_" + fHistogramNameSuffix;
  fHistConvGammaDaughtersXYRec = new TH2F(histname.Data(),"Gamma Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistNonConvGammaDaughtersXYRec_" + fHistogramNameSuffix;
  fHistNonConvGammaDaughtersXYRec = new TH2F(histname.Data(),"Gamma Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
	
  histname = "fHistNPPGammaFirstMother_" + fHistogramNameSuffix;
  fHistNPPGammaFirstMother = new TH1F(histname.Data(),"Gamma First Mother PDG Code",1201,-600.5,600.5);
  histname = "fHistNPPGammaFirstMotherXY_" + fHistogramNameSuffix;
  fHistNPPGammaFirstMotherXY = new TH2F(histname.Data(),"Gamma Mother X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistNPPGammaNDaughters_" + fHistogramNameSuffix;
  fHistNPPGammaNDaughters = new TH1F(histname.Data(),"Number of Gamma Daugthers",11,-0.5,10.5);
  histname = "fHistNPPGammaDaughters_" + fHistogramNameSuffix;
  fHistNPPGammaDaughters = new TH1F(histname.Data(),"Gamma Daugther PDG Code",1201,-600.5,600.5);
  histname = "fHistNPPGammaDaughtersXY_" + fHistogramNameSuffix;
  fHistNPPGammaDaughtersXY = new TH2F(histname.Data(),"Gamma Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
	
  histname = "fHistNPPGammaFirstMotherAcc_" + fHistogramNameSuffix;
  fHistNPPGammaFirstMotherAcc = new TH1F(histname.Data(),"Gamma First Mother PDG Code, inside EMCal acceptance",1201,-600.5,600.5);
  histname = "fHistNPPGammaFirstMotherXYAcc_" + fHistogramNameSuffix;
  fHistNPPGammaFirstMotherXYAcc = new TH2F(histname.Data(),"Gamma Mother X,Y vertex position, inside EMCal acceptance",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistNPPGammaNDaughtersAcc_" + fHistogramNameSuffix;
  fHistNPPGammaNDaughtersAcc = new TH1F(histname.Data(),"Number of Gamma Daugthers, inside EMCal acceptance",11,-0.5,10.5);
  histname = "fHistNPPGammaDaughtersAcc_" + fHistogramNameSuffix;
  fHistNPPGammaDaughtersAcc = new TH1F(histname.Data(),"Gamma Daugther PDG Code, inside EMCal acceptance",1201,-600.5,600.5);
  histname = "fHistNPPGammaDaughtersXYAcc_" + fHistogramNameSuffix;
  fHistNPPGammaDaughtersXYAcc = new TH2F(histname.Data(),"Gamma Daugther X,Y vertex position, inside EMCal acceptance",1201,-600.5,600.5,1201,-600.5,600.5);
	
  histname = "fHistNPPGammaFirstMotherRec_" + fHistogramNameSuffix;
  fHistNPPGammaFirstMotherRec = new TH1F(histname.Data(),"Reconstructed Gamma First Mother PDG Code",1201,-600.5,600.5);
  histname = "fHistNPPGammaFirstMotherXYRec_" + fHistogramNameSuffix;
  fHistNPPGammaFirstMotherXYRec = new TH2F(histname.Data(),"Gamma Mother X,Y vertex position, inside EMCal acceptance",1201,-600.5,600.5,1201,-600.5,600.5);
  histname = "fHistNPPGammaNDaughtersRec_" + fHistogramNameSuffix;
  fHistNPPGammaNDaughtersRec = new TH1F(histname.Data(),"Number of Gamma Daugthers, inside EMCal acceptance",11,-0.5,10.5);
  histname = "fHistNPPGammaDaughtersRec_" + fHistogramNameSuffix;
  fHistNPPGammaDaughtersRec = new TH1F(histname.Data(),"Gamma Daugther PDG Code",1201,-600.5,600.5);
  histname = "fHistNPPGammaDaughtersXYRec_" + fHistogramNameSuffix;
  fHistNPPGammaDaughtersXYRec = new TH2F(histname.Data(),"Gamma Daugther X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
}

void AliAnalysisEmEtMonteCarlo::FillOutputList(TList *list)
{//fill the output tlist
  //AliAnalysisEt::FillOutputList(list);

  list->Add(fHistPrimEtaEET); 
  list->Add(fHistPrimEtaPtET); 
  list->Add(fHistPrimEtaET); 
  list->Add(fHistPrimtotET); 
	
  list->Add(fHistPrimAccEtaEET); 
  list->Add(fHistPrimAccEtaPtET); 
  list->Add(fHistPrimAccEtaET); 
  list->Add(fHistPrimAcctotET); 

  list->Add(fHistPrimRecEtaEET); 
  list->Add(fHistPrimRecEtaPtET); 
  list->Add(fHistPrimRecEtaET); 
  list->Add(fHistPrimRectotET); 

  list->Add(fHistPrimRecEtaEDepETDep); 
  list->Add(fHistPrimRecEtaPtETDep); 
  list->Add(fHistPrimRecEtaETDep); 
  list->Add(fHistPrimRectotETDep); 
	
  list->Add(fHistElectronEtaEET); 
  list->Add(fHistElectronEtaPtET); 
  list->Add(fHistElectronEtaET); 
  list->Add(fHistElectronEtaE); 
  list->Add(fHistElectronEtaPt); 
  list->Add(fHistElectrontotET); 
	
  list->Add(fHistConvElectronEtaEET);  
  list->Add(fHistConvElectronEtaPtET);  
  list->Add(fHistConvElectronEtaET);  
  list->Add(fHistConvElectronEtaE);  
  list->Add(fHistConvElectronEtaPt);  
  list->Add(fHistConvElectrontotET);  
	
  list->Add(fHistScatElectronEtaEET);  
  list->Add(fHistScatElectronEtaPtET);  
  list->Add(fHistScatElectronEtaET);  
  list->Add(fHistScatElectronEtaE);  
  list->Add(fHistScatElectronEtaPt);  
  list->Add(fHistScatElectrontotET);  
	
  list->Add(fHistTotElectrontotET); 

  list->Add(fHistGammaEtaEET);  
  list->Add(fHistGammaEtaPtET);  
  list->Add(fHistGammaEtaET);  
  list->Add(fHistGammaEtaE);  
  list->Add(fHistGammaEtaPt);  
  list->Add(fHistGammatotET);  
	
  list->Add(fHistAnnihGammaEtaEET);  
  list->Add(fHistAnnihGammaEtaPtET);  
  list->Add(fHistAnnihGammaEtaET);  
  list->Add(fHistAnnihGammaEtaE);  
  list->Add(fHistAnnihGammaEtaPt);  
  list->Add(fHistAnnihGammatotET);  
	
  list->Add(fHistScatGammaEtaEET);  
  list->Add(fHistScatGammaEtaPtET);  
  list->Add(fHistScatGammaEtaET);  
  list->Add(fHistScatGammaEtaE);  
  list->Add(fHistScatGammaEtaPt);  
  list->Add(fHistScatGammatotET);  

  list->Add(fHistConvGammaEtaEET);  
  list->Add(fHistConvGammaEtaPtET);  
  list->Add(fHistConvGammaEtaET);  
  list->Add(fHistConvGammaEtaE);  
  list->Add(fHistConvGammaEtaPt);  
  list->Add(fHistConvGammatotET);  
	
  list->Add(fHistNonConvGammaEtaEET);  
  list->Add(fHistNonConvGammaEtaPtET);  
  list->Add(fHistNonConvGammaEtaET);  
  list->Add(fHistNonConvGammaEtaE);  
  list->Add(fHistNonConvGammaEtaPt);  
  list->Add(fHistNonConvGammatotET);  
		
  list->Add(fHistTotGammatotET);  

  list->Add(fHistTotEMtotET); 

  list->Add(fHistNPPElectronEtaEET); 
  list->Add(fHistNPPElectronEtaPtET); 
  list->Add(fHistNPPElectronEtaET); 
  list->Add(fHistNPPElectronEtaE); 
  list->Add(fHistNPPElectronEtaPt); 
  list->Add(fHistNPPElectrontotET); 
	
  list->Add(fHistNPPGammaEtaEET); 
  list->Add(fHistNPPGammaEtaPtET); 
  list->Add(fHistNPPGammaEtaET); 
  list->Add(fHistNPPGammaEtaE); 
  list->Add(fHistNPPGammaEtaPt); 
  list->Add(fHistNPPGammatotET); 

  list->Add(fHistTotNPPEMtotET); 

  list->Add(fHistNPPPi0GammaEtaEET); 
  list->Add(fHistNPPPi0GammaEtaPtET); 
  list->Add(fHistNPPPi0GammaEtaET); 
  list->Add(fHistNPPPi0GammaEtaE); 
  list->Add(fHistNPPPi0GammaEtaPt); 
  list->Add(fHistNPPPi0GammatotET); 
	
  list->Add(fHistElectronAccEtaEET); 
  list->Add(fHistElectronAccEtaPtET); 
  list->Add(fHistElectronAccEtaET); 
  list->Add(fHistElectronAccEtaE); 
  list->Add(fHistElectronAccEtaPt); 
  list->Add(fHistElectronAcctotET); 
	
  list->Add(fHistConvElectronAccEtaEET);  
  list->Add(fHistConvElectronAccEtaPtET);  
  list->Add(fHistConvElectronAccEtaET);  
  list->Add(fHistConvElectronAccEtaE);  
  list->Add(fHistConvElectronAccEtaPt);  
  list->Add(fHistConvElectronAcctotET);  
	
  list->Add(fHistScatElectronAccEtaEET);  
  list->Add(fHistScatElectronAccEtaPtET);  
  list->Add(fHistScatElectronAccEtaET);  
  list->Add(fHistScatElectronAccEtaE);  
  list->Add(fHistScatElectronAccEtaPt);  
  list->Add(fHistScatElectronAcctotET);  
	
  list->Add(fHistTotElectronAcctotET); 

  list->Add(fHistGammaAccEtaEET);  
  list->Add(fHistGammaAccEtaPtET);  
  list->Add(fHistGammaAccEtaET);  
  list->Add(fHistGammaAccEtaE);  
  list->Add(fHistGammaAccEtaPt);  
  list->Add(fHistGammaAcctotET);  
	
  list->Add(fHistConvGammaAccEtaEET);  
  list->Add(fHistConvGammaAccEtaPtET);  
  list->Add(fHistConvGammaAccEtaET);  
  list->Add(fHistConvGammaAccEtaE);  
  list->Add(fHistConvGammaAccEtaPt);  
  list->Add(fHistConvGammaAcctotET);  
	
  list->Add(fHistNonConvGammaAccEtaEET);  
  list->Add(fHistNonConvGammaAccEtaPtET);  
  list->Add(fHistNonConvGammaAccEtaET);  
  list->Add(fHistNonConvGammaAccEtaE);  
  list->Add(fHistNonConvGammaAccEtaPt);  
  list->Add(fHistNonConvGammaAcctotET);  
	
  list->Add(fHistAnnihGammaAccEtaEET);  
  list->Add(fHistAnnihGammaAccEtaPtET);  
  list->Add(fHistAnnihGammaAccEtaET);  
  list->Add(fHistAnnihGammaAccEtaE);  
  list->Add(fHistAnnihGammaAccEtaPt);  
  list->Add(fHistAnnihGammaAcctotET);  
	
  list->Add(fHistScatGammaAccEtaEET);  
  list->Add(fHistScatGammaAccEtaPtET);  
  list->Add(fHistScatGammaAccEtaET);  
  list->Add(fHistScatGammaAccEtaE);  
  list->Add(fHistScatGammaAccEtaPt);  
  list->Add(fHistScatGammaAcctotET);  

  list->Add(fHistTotGammaAcctotET);  

  list->Add(fHistTotEMAcctotET); 

  list->Add(fHistNPPElectronAccEtaEET); 
  list->Add(fHistNPPElectronAccEtaPtET); 
  list->Add(fHistNPPElectronAccEtaE); 
  list->Add(fHistNPPElectronAccEtaPt); 
	
  list->Add(fHistNPPGammaAccEtaEET); 
  list->Add(fHistNPPGammaAccEtaPtET); 
  list->Add(fHistNPPGammaAccEtaE); 
  list->Add(fHistNPPGammaAccEtaPt); 
	
  list->Add(fHistElectronRecEtaEET); 
  list->Add(fHistElectronRecEtaPtET); 
  list->Add(fHistElectronRecEtaET); 
  list->Add(fHistElectronRecEtaE); 
  list->Add(fHistElectronRecEtaPt); 
  list->Add(fHistElectronRectotET); 
	
  list->Add(fHistConvElectronRecEtaEET);  
  list->Add(fHistConvElectronRecEtaPtET);  
  list->Add(fHistConvElectronRecEtaET);  
  list->Add(fHistConvElectronRecEtaE);  
  list->Add(fHistConvElectronRecEtaPt);  
  list->Add(fHistConvElectronRectotET);  
	
  list->Add(fHistScatElectronRecEtaEET);  
  list->Add(fHistScatElectronRecEtaPtET);  
  list->Add(fHistScatElectronRecEtaET);  
  list->Add(fHistScatElectronRecEtaE);  
  list->Add(fHistScatElectronRecEtaPt);  
  list->Add(fHistScatElectronRectotET);  
	
  list->Add(fHistTotElectronRectotET); 

  list->Add(fHistGammaRecEtaEET);  
  list->Add(fHistGammaRecEtaPtET);  
  list->Add(fHistGammaRecEtaET);  
  list->Add(fHistGammaRecEtaE);  
  list->Add(fHistGammaRecEtaPt);  
  list->Add(fHistGammaRectotET);  
	
  list->Add(fHistAnnihGammaRecEtaEET);  
  list->Add(fHistAnnihGammaRecEtaPtET);  
  list->Add(fHistAnnihGammaRecEtaET);  
  list->Add(fHistAnnihGammaRecEtaE);  
  list->Add(fHistAnnihGammaRecEtaPt);  
  list->Add(fHistAnnihGammaRectotET);  
	
  list->Add(fHistScatGammaRecEtaEET);  
  list->Add(fHistScatGammaRecEtaPtET);  
  list->Add(fHistScatGammaRecEtaET);  
  list->Add(fHistScatGammaRecEtaE);  
  list->Add(fHistScatGammaRecEtaPt);  
  list->Add(fHistScatGammaRectotET);  

  list->Add(fHistTotGammaRectotET);  

  list->Add(fHistTotEMRectotET); 

  list->Add(fHistNPPElectronRecEtaEET); 
  list->Add(fHistNPPElectronRecEtaPtET); 
  list->Add(fHistNPPElectronRecEtaET); 
  list->Add(fHistNPPElectronRecEtaE); 
  list->Add(fHistNPPElectronRecEtaPt); 
  list->Add(fHistNPPElectronRectotET); 
	
  list->Add(fHistNPPGammaRecEtaEET); 
  list->Add(fHistNPPGammaRecEtaPtET); 
  list->Add(fHistNPPGammaRecEtaET); 
  list->Add(fHistNPPGammaRecEtaE); 
  list->Add(fHistNPPGammaRecEtaPt); 
  list->Add(fHistNPPGammaRectotET); 

  list->Add(fHistTotNPPEMRectotET); 

  list->Add(fHistNPPPi0GammaRecEtaEET); 
  list->Add(fHistNPPPi0GammaRecEtaPtET); 
  list->Add(fHistNPPPi0GammaRecEtaET); 
  list->Add(fHistNPPPi0GammaRecEtaE); 
  list->Add(fHistNPPPi0GammaRecEtaPt); 
  list->Add(fHistNPPPi0GammaRectotET); 
	
  list->Add(fHistMuonEtaEET); 
  list->Add(fHistMuonAccEtaEET); 
  list->Add(fHistMuonRecEtaEET); 
  list->Add(fHistMuonMatchEtaEET); 
	
  list->Add(fHistMuonEtaPtET); 
  list->Add(fHistMuonAccEtaPtET); 
  list->Add(fHistMuonRecEtaPtET); 
  list->Add(fHistMuonMatchEtaPtET); 
	
  list->Add(fHistMuonEtaET); 
  list->Add(fHistMuonAccEtaET); 
  list->Add(fHistMuonRecEtaET); 
  list->Add(fHistMuonMatchEtaET); 
	
  list->Add(fHistMuonEtaE); 
  list->Add(fHistMuonAccEtaE); 
  list->Add(fHistMuonRecEtaE); 
  list->Add(fHistMuonMatchEtaE); 
	
  list->Add(fHistMuonEtaPt); 
  list->Add(fHistMuonAccEtaPt); 
  list->Add(fHistMuonRecEtaPt); 
  list->Add(fHistMuonMatchEtaPt); 
	
  list->Add(fHistMuontotET); 
  list->Add(fHistMuonAcctotET); 
  list->Add(fHistMuonRectotET); 
  list->Add(fHistMuonMatchtotET); 
	
  list->Add(fHistMuonRectotETDep); 
  list->Add(fHistMuonMatchtotETDep); 
	
  list->Add(fHistMuonRecEtaEDepETDep); 
  list->Add(fHistMuonMatchEtaEDepETDep); 

  list->Add(fHistMuonRecEtaPtETDep); 
  list->Add(fHistMuonMatchEtaPtETDep); 

  list->Add(fHistMuonRecEtaETDep); 
  list->Add(fHistMuonMatchEtaETDep); 

  list->Add(fHistMuonRecResEET);
  list->Add(fHistMuonRecResPtET); 
  list->Add(fHistMuonRecResE); 
  list->Add(fHistMuonRecResPt); 
  list->Add(fHistMuonRecResEDepETDep); 
  list->Add(fHistMuonRecResPtETDep); 
	
  list->Add(fHistPionEtaEET); 
  list->Add(fHistPionAccEtaEET); 
  list->Add(fHistPionRecEtaEET); 
  list->Add(fHistPionMatchEtaEET); 
	
  list->Add(fHistPionEtaPtET); 
  list->Add(fHistPionAccEtaPtET); 
  list->Add(fHistPionRecEtaPtET); 
  list->Add(fHistPionMatchEtaPtET); 
	
  list->Add(fHistPionEtaET); 
  list->Add(fHistPionAccEtaET); 
  list->Add(fHistPionRecEtaET); 
  list->Add(fHistPionMatchEtaET); 
	
  list->Add(fHistPionEtaE); 
  list->Add(fHistPionAccEtaE); 
  list->Add(fHistPionRecEtaE); 
  list->Add(fHistPionMatchEtaE); 
	
  list->Add(fHistPionEtaPt); 
  list->Add(fHistPionAccEtaPt); 
  list->Add(fHistPionRecEtaPt); 
  list->Add(fHistPionMatchEtaPt); 
	
  list->Add(fHistPiontotET); 
  list->Add(fHistPionAcctotET); 
  list->Add(fHistPionRectotET); 
  list->Add(fHistPionMatchtotET); 
	
  list->Add(fHistPionRectotETDep); 
  list->Add(fHistPionMatchtotETDep); 
	
  list->Add(fHistPionRecEtaEDepETDep); 
  list->Add(fHistPionMatchEtaEDepETDep); 
	
  list->Add(fHistPionRecEtaPtETDep); 
  list->Add(fHistPionMatchEtaPtETDep); 
	
  list->Add(fHistPionRecEtaETDep); 
  list->Add(fHistPionMatchEtaETDep); 
	
  list->Add(fHistPionRecResEET);
  list->Add(fHistPionRecResPtET); 
  list->Add(fHistPionRecResE); 
  list->Add(fHistPionRecResPt); 
  list->Add(fHistPionRecResEDepETDep); 
  list->Add(fHistPionRecResPtETDep); 
	
  list->Add(fHistKaonEtaEET); 
  list->Add(fHistKaonAccEtaEET); 
  list->Add(fHistKaonRecEtaEET); 
  list->Add(fHistKaonMatchEtaEET); 
	
  list->Add(fHistKaonEtaPtET); 
  list->Add(fHistKaonAccEtaPtET); 
  list->Add(fHistKaonRecEtaPtET); 
  list->Add(fHistKaonMatchEtaPtET); 
	
  list->Add(fHistKaonEtaET); 
  list->Add(fHistKaonAccEtaET); 
  list->Add(fHistKaonRecEtaET); 
  list->Add(fHistKaonMatchEtaET); 
	
  list->Add(fHistKaonEtaE); 
  list->Add(fHistKaonAccEtaE); 
  list->Add(fHistKaonRecEtaE); 
  list->Add(fHistKaonMatchEtaE); 
	
  list->Add(fHistKaonEtaPt); 
  list->Add(fHistKaonAccEtaPt); 
  list->Add(fHistKaonRecEtaPt); 
  list->Add(fHistKaonMatchEtaPt); 
	
  list->Add(fHistKaontotET); 
  list->Add(fHistKaonAcctotET); 
  list->Add(fHistKaonRectotET); 
  list->Add(fHistKaonMatchtotET); 
	
  list->Add(fHistKaonRectotETDep); 
  list->Add(fHistKaonMatchtotETDep); 
	
  list->Add(fHistKaonRecEtaEDepETDep); 
  list->Add(fHistKaonMatchEtaEDepETDep); 
	
  list->Add(fHistKaonRecEtaPtETDep); 
  list->Add(fHistKaonMatchEtaPtETDep); 
	
  list->Add(fHistKaonRecEtaETDep); 
  list->Add(fHistKaonMatchEtaETDep); 
	
  list->Add(fHistKaonRecResEET);
  list->Add(fHistKaonRecResPtET); 
  list->Add(fHistKaonRecResE); 
  list->Add(fHistKaonRecResPt); 
  list->Add(fHistKaonRecResEDepETDep); 
  list->Add(fHistKaonRecResPtETDep); 
	
  list->Add(fHistProtonEtaEET); 
  list->Add(fHistProtonAccEtaEET); 
  list->Add(fHistProtonRecEtaEET); 
  list->Add(fHistProtonMatchEtaEET); 
	
  list->Add(fHistProtonEtaPtET); 
  list->Add(fHistProtonAccEtaPtET); 
  list->Add(fHistProtonRecEtaPtET); 
  list->Add(fHistProtonMatchEtaPtET); 
	
  list->Add(fHistProtonEtaET); 
  list->Add(fHistProtonAccEtaET); 
  list->Add(fHistProtonRecEtaET); 
  list->Add(fHistProtonMatchEtaET); 
	
  list->Add(fHistProtonEtaE); 
  list->Add(fHistProtonAccEtaE); 
  list->Add(fHistProtonRecEtaE); 
  list->Add(fHistProtonMatchEtaE); 
	
  list->Add(fHistProtonEtaPt); 
  list->Add(fHistProtonAccEtaPt); 
  list->Add(fHistProtonRecEtaPt); 
  list->Add(fHistProtonMatchEtaPt); 

  list->Add(fHistProtontotET); 
  list->Add(fHistProtonAcctotET); 
  list->Add(fHistProtonRectotET); 
  list->Add(fHistProtonMatchtotET); 
	
  list->Add(fHistProtonRectotETDep); 
  list->Add(fHistProtonMatchtotETDep); 
	
  list->Add(fHistProtonRecEtaEDepETDep); 
  list->Add(fHistProtonMatchEtaEDepETDep); 
	
  list->Add(fHistProtonRecEtaPtETDep); 
  list->Add(fHistProtonMatchEtaPtETDep); 
	
  list->Add(fHistProtonRecEtaETDep); 
  list->Add(fHistProtonMatchEtaETDep); 
	
  list->Add(fHistProtonRecResEET);
  list->Add(fHistProtonRecResPtET); 
  list->Add(fHistProtonRecResE); 
  list->Add(fHistProtonRecResPt); 
  list->Add(fHistProtonRecResEDepETDep); 
  list->Add(fHistProtonRecResPtETDep); 

  list->Add(fHistTotChargedtotET); 
  list->Add(fHistTotChargedAcctotET); 
  list->Add(fHistTotChargedRectotET); 
  list->Add(fHistTotChargedMatchtotET); 
	
  list->Add(fHistTotChargedRectotETDep); 
  list->Add(fHistTotChargedMatchtotETDep); 
	
  list->Add(fHistNeutronEtaEET); 
  list->Add(fHistNeutronAccEtaEET); 
  list->Add(fHistNeutronRecEtaEET); 
	
  list->Add(fHistNeutronEtaPtET); 
  list->Add(fHistNeutronAccEtaPtET); 
  list->Add(fHistNeutronRecEtaPtET); 
	
  list->Add(fHistNeutronEtaET); 
  list->Add(fHistNeutronAccEtaET); 
  list->Add(fHistNeutronRecEtaET); 
	
  list->Add(fHistNeutronEtaE); 
  list->Add(fHistNeutronAccEtaE); 
  list->Add(fHistNeutronRecEtaE); 
	
  list->Add(fHistNeutronEtaPt); 
  list->Add(fHistNeutronAccEtaPt); 
  list->Add(fHistNeutronRecEtaPt); 
	
  list->Add(fHistNeutrontotET); 
  list->Add(fHistNeutronAcctotET); 
  list->Add(fHistNeutronRectotET); 
  list->Add(fHistNeutronRectotETDep); 
	
  list->Add(fHistNeutronRecEtaEDepETDep); 
  list->Add(fHistNeutronRecEtaETDep); 
		
  list->Add(fHistNeutronRecEtaPtETDep); 
	
  list->Add(fHistK0EtaEET); 
  list->Add(fHistK0RecEtaEET); 
	
  list->Add(fHistK0EtaPtET); 
  list->Add(fHistK0RecEtaPtET); 
	
  list->Add(fHistK0EtaET); 
  list->Add(fHistK0RecEtaET); 
	
  list->Add(fHistK0EtaE); 
  list->Add(fHistK0RecEtaE); 
	
  list->Add(fHistK0EtaPt); 
  list->Add(fHistK0RecEtaPt); 
	
  list->Add(fHistK0totET); 
  list->Add(fHistK0RectotET); 
  list->Add(fHistK0RectotETDep); 
	
  list->Add(fHistK0RecEtaEDepETDep); 
  list->Add(fHistK0RecEtaETDep); 
	
  list->Add(fHistK0RecEtaPtETDep); 
	
  list->Add(fHistLambdaEtaEET); 
  list->Add(fHistLambdaRecEtaEET); 
	
  list->Add(fHistLambdaEtaPtET); 
  list->Add(fHistLambdaRecEtaPtET); 
	
  list->Add(fHistLambdaEtaET); 
  list->Add(fHistLambdaRecEtaET); 
	
  list->Add(fHistLambdaEtaE); 
  list->Add(fHistLambdaRecEtaE); 
	
  list->Add(fHistLambdaEtaPt); 
  list->Add(fHistLambdaRecEtaPt); 
	
  list->Add(fHistLambdatotET); 
  list->Add(fHistLambdaRectotET); 
  list->Add(fHistLambdaRectotETDep); 
	
  list->Add(fHistLambdaRecEtaEDepETDep); 
  list->Add(fHistLambdaRecEtaETDep); 
	
  list->Add(fHistLambdaRecEtaPtETDep); 
	
  list->Add(fHistTotNeutraltotET); 
  list->Add(fHistTotNeutralRectotET); 
  list->Add(fHistTotNeutralRectotETDep); 

  list->Add(fHistTotaltotET); 
  list->Add(fHistTotalAcctotET); 
  list->Add(fHistTotalRectotET); 
  list->Add(fHistTotalRectotETDep); 
	
  list->Add(fHistElectronFirstMother); 
  list->Add(fHistElectronFirstMotherXY); 
  list->Add(fHistElectronNDaughters); 
  list->Add(fHistElectronDaughters); 
  list->Add(fHistElectronDaughtersXY); 
	
  list->Add(fHistElectronFirstMotherAcc);  
  list->Add(fHistElectronFirstMotherXYAcc);  
  list->Add(fHistElectronNDaughtersAcc); 
  list->Add(fHistElectronDaughtersAcc); 
  list->Add(fHistElectronDaughtersXYAcc); 
	
  list->Add(fHistElectronFirstMotherRec);  
  list->Add(fHistElectronFirstMotherXYRec);  
  list->Add(fHistElectronNDaughtersRec); 
  list->Add(fHistElectronDaughtersRec); 
  list->Add(fHistElectronDaughtersXYRec); 
	
  list->Add(fHistNPPElectronFirstMother); 
  list->Add(fHistNPPElectronFirstMotherXY); 
  list->Add(fHistNPPElectronNDaughters); 
  list->Add(fHistNPPElectronDaughters); 
  list->Add(fHistNPPElectronDaughtersXY); 
	
  list->Add(fHistNPPElectronFirstMotherAcc);  
  list->Add(fHistNPPElectronFirstMotherXYAcc);  
  list->Add(fHistNPPElectronNDaughtersAcc); 
  list->Add(fHistNPPElectronDaughtersAcc); 
  list->Add(fHistNPPElectronDaughtersXYAcc); 
	
  list->Add(fHistNPPElectronFirstMotherRec);  
  list->Add(fHistNPPElectronFirstMotherXYRec);  
  list->Add(fHistNPPElectronNDaughtersRec); 
  list->Add(fHistNPPElectronDaughtersRec); 
  list->Add(fHistNPPElectronDaughtersXYRec); 
	
  list->Add(fHistGammaFirstMother); 
  list->Add(fHistGammaFirstMotherXY); 
  list->Add(fHistGammaNDaughters); 
  list->Add(fHistGammaDaughters); 
  list->Add(fHistGammaDaughtersXY); 
  list->Add(fHistConvGammaDaughtersXY); 
  list->Add(fHistNonConvGammaDaughtersXY); 
	
  list->Add(fHistGammaFirstMotherAcc);  
  list->Add(fHistGammaFirstMotherXYAcc);  
  list->Add(fHistGammaNDaughtersAcc); 
  list->Add(fHistGammaDaughtersAcc); 
  list->Add(fHistGammaDaughtersXYAcc); 
  list->Add(fHistConvGammaDaughtersXYAcc); 
  list->Add(fHistNonConvGammaDaughtersXYAcc); 
	
  list->Add(fHistGammaFirstMotherRec);  
  list->Add(fHistGammaFirstMotherXYRec);  
  list->Add(fHistGammaNDaughtersRec); 
  list->Add(fHistGammaDaughtersRec); 
  list->Add(fHistGammaDaughtersXYRec); 
  list->Add(fHistConvGammaDaughtersXYRec); 
  list->Add(fHistNonConvGammaDaughtersXYRec); 
	
  list->Add(fHistNPPGammaFirstMother); 
  list->Add(fHistNPPGammaFirstMotherXY); 
  list->Add(fHistNPPGammaNDaughters); 
  list->Add(fHistNPPGammaDaughters); 
  list->Add(fHistNPPGammaDaughtersXY); 
	
  list->Add(fHistNPPGammaFirstMotherAcc);  
  list->Add(fHistNPPGammaFirstMotherXYAcc);  
  list->Add(fHistNPPGammaNDaughtersAcc); 
  list->Add(fHistNPPGammaDaughtersAcc); 
  list->Add(fHistNPPGammaDaughtersXYAcc); 
	
  list->Add(fHistNPPGammaFirstMotherRec);  
  list->Add(fHistNPPGammaFirstMotherXYRec);  
  list->Add(fHistNPPGammaNDaughtersRec); 
  list->Add(fHistNPPGammaDaughtersRec); 
  list->Add(fHistNPPGammaDaughtersXYRec); 
	
  list->Add(fHistAllERecEMC);	
  list->Add(fHistAllPtRecPtMC);
  list->Add(fHistElectronERecEMC);	
  list->Add(fHistGammaERecEMC);
	
  list->Add(fHistChargedRes);
  list->Add(fHistChargedRes2);
  list->Add(fHistChargedRes3);
  list->Add(fHistNeutralRes);
  list->Add(fHistElectronRes);
  list->Add(fHistGammaRes);
	
  list->Add(fHistIsInAcc);
}

//________________________________________________________________________
Bool_t AliAnalysisEmEtMonteCarlo::TrackHitsCalo(AliExternalTrackParam* extParam)
{//Does the track hit the calorimeter?
  TVector3 pos(0,0,0);
	
  if (extParam)
    {
      if (GetTrackProjection(extParam,pos))
	{
	  Bool_t inAcc = fGeoUt->IsInEMCAL(pos.X(),pos.Y(),pos.Z());
			
	  //if (inAcc)
	  //	fHistIsInAcc->Fill(pos.X(),pos.Y());
			
	  return inAcc;
	}
    }
  return kFALSE;
}

//________________________________________________________________________
//project to a EMCal radius
Bool_t AliAnalysisEmEtMonteCarlo::GetTrackProjection(AliExternalTrackParam *trackParam, TVector3 &trackPos)
{//Get the track projection
  Bool_t proj = kFALSE;
  Double_t emcalR = fGeoUt->GetEMCGeometry()->GetIPDistance();
	
  if (trackParam) //it is constructed from TParticle
    {
      Double_t trkPos[3] = {0};
		
      //Assume the track is a pion with mass 0.139GeV/c^2
      //Extrapolation step is 1cm
      if(!AliTrackerBase::PropagateTrackToBxByBz(trackParam, emcalR, 0.139, 1, kTRUE, 0.8) ) return proj;
		
      trackParam->GetXYZ(trkPos);
		
      trackPos.SetXYZ(trkPos[0],trkPos[1],trkPos[2]);
		
      proj = kTRUE;               
    }
	
  return proj;
}

//________________________________________________________________________
//project to a cluster position
Bool_t AliAnalysisEmEtMonteCarlo::GetTrackProjection(AliEMCALTrack* emcTrack, TVector3 &trackPos, TVector3 clusPos)
{//get the track projection
  Bool_t proj = kFALSE;
	
  if (emcTrack)
    {	
      Double_t trkPos[3] = {0};
		
      emcTrack->PropagateToGlobal(clusPos.X(),clusPos.Y(),clusPos.Z(),0.,0.);
      emcTrack->GetXYZ(trkPos);
		
      trackPos.SetXYZ(trkPos[0],trkPos[1],trkPos[2]);
		
      proj = kTRUE;
    }
	
  return proj;
}

//________________________________________________________________________
Bool_t AliAnalysisEmEtMonteCarlo::IsInAcceptance(TParticle *part, TParticlePDG *pdg, AliExternalTrackParam* extParam)
{//is the track in the acceptance of the emcal?
  if ((part) && (pdg))
    {
      if (TMath::Abs(pdg->Charge() - fCuts->GetMonteCarloNeutralParticle()) <1e-3 )
	{
	  if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiCutAccMax && part->Phi() > fPhiCutAccMin) 
	    return kTRUE;
	}
      else
	{
	  return TrackHitsCalo(extParam);
	}
    }

  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisEmEtMonteCarlo::IsInAcceptance(AliMCParticle *part)
{// is the track in the acceptance of the emcal?
  if (part)
    {
      for (int i=0;i<part->GetNumberOfTrackReferences();i++)
	{
	  AliTrackReference* aliTrkRef = part->GetTrackReference(i);
			
	  if (aliTrkRef)
	    {
	      //if (aliTrkRef->DetectorId() == AliTrackReference::kEMCAL)
	      //	return kTRUE;
	      if ( (aliTrkRef->DetectorId() == AliTrackReference::kEMCAL) || (fGeoUt->IsInEMCAL(aliTrkRef->X(),aliTrkRef->Y(),aliTrkRef->Z())) )
		return kTRUE;
	    }
	}
    }
	
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisEmEtMonteCarlo::IsPrimary(AliStack *stack, Int_t iPart, TParticlePDG *pdg, Int_t iPartMom, TParticlePDG *pdgMom)
{//Is the track a primary track?
  if (stack->IsPhysicalPrimary(iPart))
    {
      return kTRUE;
    }
  else if (pdg)
    {
      if (((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode) || (pdg->PdgCode() == fgGammaCode)) && 
	  ((IsMotherPrimaryGamma(stack,iPartMom,pdgMom)) || (IsMotherPrimaryElectron(stack,iPartMom,pdgMom)))  )
	{
	  return kTRUE;
	}
    }

  return kFALSE;	
}

//________________________________________________________________________
Bool_t AliAnalysisEmEtMonteCarlo::IsMotherPrimaryGamma(AliStack *stack, Int_t iPartMom, TParticlePDG *pdgMom)
{//Is the mother a primary gamma?
  Int_t nStackTracks = stack->GetNtrack();

  if (pdgMom)
    {
      if ((pdgMom->PdgCode() == fgGammaCode) && (iPartMom>=0) && (iPartMom < nStackTracks))
	{
	  if (stack->IsPhysicalPrimary(iPartMom))
	    return kTRUE;
	}
    }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisEmEtMonteCarlo::IsMotherPrimaryElectron(AliStack *stack, Int_t iPartMom, TParticlePDG *pdgMom)
{//is the mother a primary electron?
  Int_t nStackTracks = stack->GetNtrack();

  if (pdgMom)
    {
      if ((pdgMom->PdgCode() == fgEPlusCode || pdgMom->PdgCode() == fgEMinusCode) && (iPartMom>=0) && (iPartMom < nStackTracks))
	{
	  if (stack->IsPhysicalPrimary(iPartMom))
	    return kTRUE;
	}
    }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisEmEtMonteCarlo::IsGammaConversion(AliStack *stack, TParticle *part, TParticlePDG *pdg)
{//is the gamma converted?
  if ((part) && (pdg))
    {
      if(pdg->PdgCode() == fgGammaCode)
	{
	  Int_t nStackTracks = stack->GetNtrack();
	  Int_t iPartDaughter = part->GetLastDaughter();
			
	  if ((iPartDaughter>=0) && (iPartDaughter < nStackTracks))
	    {
	      TParticle *partDaughter = stack->Particle(iPartDaughter);
	      if (partDaughter)
		{
		  TParticlePDG *pdgDaugther = partDaughter->GetPDG(0);
		  if (pdgDaugther) 
		    {
		      if ( ((pdgDaugther->PdgCode() == fgEPlusCode) || (pdgDaugther->PdgCode() == fgEMinusCode)) && (!fGeoUt->IsInEMCAL(partDaughter->Vx(),partDaughter->Vy(),partDaughter->Vz())) )
			{
			  //Double_t emcalR = fGeoUt->GetEMCGeometry()->GetIPDistance();
			  //Double_t decayR = sqrt(pow(partDaughter->Vx(),2)+pow(partDaughter->Vy(),2));
							
			  //if (decayR<emcalR)						
			  return kTRUE;
			}
		    }
		}
	    }
	}
    }
  return kFALSE;			
}
	
//________________________________________________________________________
AliExternalTrackParam* AliAnalysisEmEtMonteCarlo::CreateExternalTrackParam(TParticle *part)
{//create external track param
  // Calculate the AliExternalTrackParam content
  Double_t xref;
  Double_t alpha;
  Double_t param[5];
  Double_t covar[15];
	
  // Calculate alpha: the rotation angle of the corresponding local system (TPC sector)
  alpha = part->Phi()*180./TMath::Pi();
  if (alpha<0) alpha+= 360.;
  if (alpha>360) alpha -= 360.;
	
  Int_t sector = (Int_t)(alpha/20.);
  alpha = 10. + 20.*sector;
  alpha /= 180;
  alpha *= TMath::Pi();
	
  // Covariance matrix: no errors, the parameters are exact
  for (int i=0; i<15; i++) covar[i]=0.;
	
  // Get the vertex of origin and the momentum
  TVector3 ver(part->Vx(),part->Vy(),part->Vz());
  TVector3 mom(part->Px(),part->Py(),part->Pz());
	
  // Rotate to the local coordinate system (TPC sector)
  ver.RotateZ(-alpha);
  mom.RotateZ(-alpha);
	
  // X of the referense plane
  xref = ver.X();
	
  Double_t charge;
  if (part->GetPDG(0))
    charge = part->GetPDG(0)->Charge();
  else
    return 0;
	
  if (mom.Pt()>0)
    {
      param[0] = ver.Y();
      param[1] = ver.Z();
      param[2] = TMath::Sin(mom.Phi());
      param[3] = mom.Pz()/mom.Pt();
      param[4] = TMath::Sign(1/mom.Pt(),charge);
    }
  else
    return 0;
	
  // Set AliExternalTrackParam
  AliExternalTrackParam* extTrkParam = new AliExternalTrackParam(xref, alpha, param, covar);
	
  return extTrkParam;
}

//________________________________________________________________________
Double_t AliAnalysisEmEtMonteCarlo::CalcET(TParticle *part, TParticlePDG *pdg)
{//Calculate Et
  //***************
    // calculate E_T
    //***************
    Double_t particleMassPart = 0; //The mass part in the Et calculation for this particle
  Double_t protonMass = fgProtonMass;

  if (pdg)
    {
      if (
	  TMath::Abs(pdg->PdgCode()) == fgProtonCode ||
	  TMath::Abs(pdg->PdgCode()) == fgNeutronCode ||
	  TMath::Abs(pdg->PdgCode()) == fgLambdaCode ||
	  TMath::Abs(pdg->PdgCode()) == fgXiCode ||
	  TMath::Abs(pdg->PdgCode()) == fgXi0Code ||
	  TMath::Abs(pdg->PdgCode()) == fgOmegaCode
	  )
	{
	  if (pdg->PdgCode() > 0) { particleMassPart = - protonMass;}
	  if (pdg->PdgCode() < 0) { particleMassPart = protonMass;}
	}
      Double_t et = part->Energy() * TMath::Sin(part->Theta()) + particleMassPart;
      return et;
    }
  else
    return -1.;
}
	
//________________________________________________________________________
Double_t AliAnalysisEmEtMonteCarlo::CalcETDep(Double_t caloE, TParticle *part, TParticlePDG *pdg)
{//calculate et dependence
  //***************
    // calculate E_T
    //***************
    Double_t particleMassPart = 0; //The mass part in the Et calculation for this particle
  Double_t protonMass = fgProtonMass;
	
  if (pdg)
    {
      if (
	  TMath::Abs(pdg->PdgCode()) == fgProtonCode ||
	  TMath::Abs(pdg->PdgCode()) == fgNeutronCode ||
	  TMath::Abs(pdg->PdgCode()) == fgLambdaCode ||
	  TMath::Abs(pdg->PdgCode()) == fgXiCode ||
	  TMath::Abs(pdg->PdgCode()) == fgXi0Code ||
	  TMath::Abs(pdg->PdgCode()) == fgOmegaCode
	  )
	{
	  if (pdg->PdgCode() > 0) { particleMassPart = - protonMass;}
	  if (pdg->PdgCode() < 0) { particleMassPart = protonMass;}
	}
      Double_t et = caloE * TMath::Sin(part->Theta()) + particleMassPart;
      return et;
    }
  else
    return -1.;
}


	
