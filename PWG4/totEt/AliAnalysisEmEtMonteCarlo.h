#ifndef ALIANALYSISEMETMONTECARLO_H
#define ALIANALYSISEMETMONTECARLO_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis
//  - MC output
//
//*-- Author: Marcelo G. Munhoz (USP)
//_________________________________________________________________________

#include "AliAnalysisEtMonteCarlo.h"
class TParticle;
class TParticlePDG;
class AliMCParticle;
class AliESDtrack;
class AliEMCALTrack;
class TVector3;
class AliEMCALGeometry;
class AliExternalTrackParam;
class AliStack;

class AliAnalysisEmEtMonteCarlo : public AliAnalysisEtMonteCarlo
{

public:
   
  AliAnalysisEmEtMonteCarlo();
  virtual ~AliAnalysisEmEtMonteCarlo();

    virtual Int_t AnalyseEvent(AliVEvent* event);
	virtual Int_t AnalyseEvent(AliVEvent* event, AliVEvent* event2);

    virtual void Init();
    virtual void ResetEventValues();
    virtual void CreateHistograms();
    virtual void FillOutputList(TList* list);

protected:

	virtual Bool_t IsPrimary(AliStack *stack, Int_t part, TParticlePDG *pdg, Int_t partMom, TParticlePDG *pdgMom);
	virtual Bool_t IsMotherPrimaryGamma(AliStack *stack, Int_t iPartMom, TParticlePDG *pdgMom);
	virtual Bool_t IsMotherPrimaryElectron(AliStack *stack, Int_t iPartMom, TParticlePDG *pdgMom);
	virtual Bool_t IsGammaConversion(AliStack *stack, TParticle *part, TParticlePDG *pdg);
	virtual Bool_t IsInAcceptance(TParticle *part=0, TParticlePDG *pdg=0, AliExternalTrackParam* extParam=0);
	virtual Bool_t IsInAcceptance(AliMCParticle *part=0);
	
	virtual Bool_t TrackHitsCalo(AliExternalTrackParam *extParam);
	
	virtual Bool_t GetTrackProjection(AliExternalTrackParam *trackParam, TVector3 &trackPos); // project to a radius
	virtual Bool_t GetTrackProjection(AliEMCALTrack* emcTrack, TVector3 &trackPos, TVector3 clusPos); // project to a point

	AliExternalTrackParam* CreateExternalTrackParam(TParticle *part);
	
	virtual Double_t CalcET(TParticle *part, TParticlePDG *pdg);
	virtual Double_t CalcETDep(Double_t caloE, TParticle *part, TParticlePDG *pdg);
	
protected:

	//Int_t fNcoll; // Ncoll, for Hijing; 1 otherwise
	//Int_t fNpart; // Ncoll, for Hijing; 2 otherwise

	//Double_t fImpactParameter, fResCut; // b(fm), for Hijing; 0 otherwise
	Double_t fResCut; // b(fm), for Hijing; 0 otherwise
	Double_t fPrimtotET, fPrimAcctotET, fPrimRectotET, fPrimRectotETDep;//Marcelo please add comment

	Double_t fElectrontotET, fElectronAcctotET, fElectronRectotET;//Marcelo please add comment
	Double_t fConvElectrontotET, fConvElectronAcctotET, fConvElectronRectotET, fScatElectrontotET, fScatElectronAcctotET, fScatElectronRectotET;//Marcelo please add comment
	Double_t fTotElectrontotET, fTotElectronAcctotET, fTotElectronRectotET;//Marcelo please add comment

	Double_t fGammatotET, fGammaAcctotET, fGammaRectotET;//Marcelo please add comment
	Double_t fAnnihGammatotET, fAnnihGammaAcctotET, fAnnihGammaRectotET, fScatGammatotET, fScatGammaAcctotET, fScatGammaRectotET;//Marcelo please add comment
	Double_t fTotGammatotET, fTotGammaAcctotET, fTotGammaRectotET;//Marcelo please add comment
	Double_t fConvGammatotET, fNonConvGammatotET, fConvGammaAcctotET, fNonConvGammaAcctotET, fNPPPi0GammatotET, fNPPPi0GammaRectotET;//Marcelo please add comment

	Double_t fTotEMtotET, fTotEMAcctotET, fTotEMRectotET;//Marcelo please add comment

	Double_t fNPPElectrontotET, fNPPElectronRectotET, fNPPGammatotET, fNPPGammaRectotET;//Marcelo please add comment
	Double_t fTotNPPEMtotET, fTotNPPEMRectotET;//Marcelo please add comment

	Double_t fMuontotET, fPiontotET, fKaontotET, fProtontotET;//Marcelo please add comment
	Double_t fMuonAcctotET, fPionAcctotET, fKaonAcctotET, fProtonAcctotET;//Marcelo please add comment
	Double_t fMuonRectotET, fMuonRectotETDep, fPionRectotET, fPionRectotETDep, fKaonRectotET, fKaonRectotETDep, fProtonRectotET, fProtonRectotETDep;//Marcelo please add comment
	Double_t fMuonMatchtotET, fMuonMatchtotETDep, fPionMatchtotET, fPionMatchtotETDep, fKaonMatchtotET, fKaonMatchtotETDep, fProtonMatchtotET, fProtonMatchtotETDep;//Marcelo please add comment
	Double_t fTotChargedtotET, fTotChargedAcctotET, fTotChargedRectotET, fTotChargedRectotETDep, fTotChargedMatchtotET, fTotChargedMatchtotETDep;//Marcelo please add comment

	Double_t fNeutrontotET, fNeutronAcctotET, fNeutronRectotET, fNeutronRectotETDep;//Marcelo please add comment
	Double_t fK0totET, fK0RectotET, fK0RectotETDep, fLambdatotET, fLambdaRectotET, fLambdaRectotETDep;//Marcelo please add comment
	Double_t fTotNeutraltotET, fTotNeutralRectotET, fTotNeutralRectotETDep;//Marcelo please add comment

	Double_t fTotaltotET, fTotalAcctotET, fTotalRectotET, fTotalRectotETDep;//Marcelo please add comment
	
	AliEMCALGeometry *fGeoUt;//Marcelo please add comment

	// *******************
	// primaries ET
	// *******************
	TH2F *fHistPrimEtaEET;//Marcelo please add comment 
	TH2F *fHistPrimEtaPtET;//Marcelo please add comment 
	TH2F *fHistPrimEtaET;//Marcelo please add comment
	TH1F *fHistPrimtotET;//Marcelo please add comment
	
	TH2F *fHistPrimAccEtaEET;//Marcelo please add comment 
	TH2F *fHistPrimAccEtaPtET;//Marcelo please add comment 
	TH2F *fHistPrimAccEtaET;//Marcelo please add comment 
	TH1F *fHistPrimAcctotET;//Marcelo please add comment
	
	TH2F *fHistPrimRecEtaEET;//Marcelo please add comment 
	TH2F *fHistPrimRecEtaPtET;//Marcelo please add comment 
	TH2F *fHistPrimRecEtaET;//Marcelo please add comment 
	TH1F *fHistPrimRectotET;//Marcelo please add comment

	TH2F *fHistPrimRecEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistPrimRecEtaPtETDep;//Marcelo please add comment 
	TH2F *fHistPrimRecEtaETDep;//Marcelo please add comment 
	TH1F *fHistPrimRectotETDep;//Marcelo please add comment
	
	// *******************
	// electron ET
	// *******************
	TH2F *fHistElectronEtaEET;//Marcelo please add comment 
	TH2F *fHistElectronEtaPtET;//Marcelo please add comment 
	TH2F *fHistElectronEtaET;//Marcelo please add comment 
	TH2F *fHistElectronEtaE;//Marcelo please add comment 
	TH2F *fHistElectronEtaPt;//Marcelo please add comment 
	TH1F *fHistElectrontotET;//Marcelo please add comment 

	//TH2F *fHistConvElectronEtaEET;//Marcelo please add comment  
	TH2F *fHistConvElectronEtaPtET;//Marcelo please add comment  
	TH2F *fHistConvElectronEtaET;//Marcelo please add comment  
	//TH2F *fHistConvElectronEtaE;//Marcelo please add comment  
	TH2F *fHistConvElectronEtaPt;//Marcelo please add comment  
	TH1F *fHistConvElectrontotET;//Marcelo please add comment  

	TH2F *fHistScatElectronEtaEET;//Marcelo please add comment  
	TH2F *fHistScatElectronEtaPtET;//Marcelo please add comment  
	TH2F *fHistScatElectronEtaET;//Marcelo please add comment  
	TH2F *fHistScatElectronEtaE;//Marcelo please add comment  
	TH2F *fHistScatElectronEtaPt;//Marcelo please add comment  
	TH1F *fHistScatElectrontotET;//Marcelo please add comment  
	
	// *******************
	// total electron ET
	// *******************
	TH1F *fHistTotElectrontotET;//Marcelo please add comment
	
	// *******************
	// gamma ET
	// *******************
	TH2F *fHistGammaEtaEET;//Marcelo please add comment  
	TH2F *fHistGammaEtaPtET;//Marcelo please add comment  
	TH2F *fHistGammaEtaET;//Marcelo please add comment  
	TH2F *fHistGammaEtaE;//Marcelo please add comment  
	TH2F *fHistGammaEtaPt;//Marcelo please add comment  
	TH1F *fHistGammatotET;//Marcelo please add comment  
	
	TH2F *fHistAnnihGammaEtaEET;//Marcelo please add comment  
	TH2F *fHistAnnihGammaEtaPtET;//Marcelo please add comment  
	TH2F *fHistAnnihGammaEtaET;//Marcelo please add comment  
	TH2F *fHistAnnihGammaEtaE;//Marcelo please add comment  
	TH2F *fHistAnnihGammaEtaPt;//Marcelo please add comment  
	TH1F *fHistAnnihGammatotET;//Marcelo please add comment  

	TH2F *fHistScatGammaEtaEET;//Marcelo please add comment  
	TH2F *fHistScatGammaEtaPtET;//Marcelo please add comment  
	TH2F *fHistScatGammaEtaET;//Marcelo please add comment  
	TH2F *fHistScatGammaEtaE;//Marcelo please add comment  
	TH2F *fHistScatGammaEtaPt;//Marcelo please add comment  
	TH1F *fHistScatGammatotET;//Marcelo please add comment  

	TH2F *fHistConvGammaEtaEET;//Marcelo please add comment  
	TH2F *fHistConvGammaEtaPtET;//Marcelo please add comment  
	TH2F *fHistConvGammaEtaET;//Marcelo please add comment  
	TH2F *fHistConvGammaEtaE;//Marcelo please add comment  
	TH2F *fHistConvGammaEtaPt;//Marcelo please add comment  
	TH1F *fHistConvGammatotET;//Marcelo please add comment  
	
	TH2F *fHistNonConvGammaEtaEET;//Marcelo please add comment  
	TH2F *fHistNonConvGammaEtaPtET;//Marcelo please add comment  
	TH2F *fHistNonConvGammaEtaET;//Marcelo please add comment  
	TH2F *fHistNonConvGammaEtaE;//Marcelo please add comment  
	TH2F *fHistNonConvGammaEtaPt;//Marcelo please add comment  
	TH1F *fHistNonConvGammatotET;//Marcelo please add comment  
	
	// *******************
	// total gamma ET
	// *******************
	TH1F *fHistTotGammatotET;//Marcelo please add comment

	// *******************
	// total electromagnetic ET
	// *******************
	TH1F *fHistTotEMtotET;//Marcelo please add comment

	// non-primary electromagnetic ET
	TH2F *fHistNPPElectronEtaEET;//Marcelo please add comment 
	TH2F *fHistNPPElectronEtaPtET;//Marcelo please add comment 
	TH2F *fHistNPPElectronEtaET;//Marcelo please add comment 
	TH2F *fHistNPPElectronEtaE;//Marcelo please add comment 
	TH2F *fHistNPPElectronEtaPt;//Marcelo please add comment 
	TH1F *fHistNPPElectrontotET;//Marcelo please add comment 

	TH2F *fHistNPPGammaEtaEET;//Marcelo please add comment 
	TH2F *fHistNPPGammaEtaPtET;//Marcelo please add comment 
	TH2F *fHistNPPGammaEtaET;//Marcelo please add comment 
	TH2F *fHistNPPGammaEtaE;//Marcelo please add comment 
	TH2F *fHistNPPGammaEtaPt;//Marcelo please add comment 
	TH1F *fHistNPPGammatotET;//Marcelo please add comment 

	TH1F *fHistTotNPPEMtotET;//Marcelo please add comment

	TH2F *fHistNPPPi0GammaEtaEET;//Marcelo please add comment 
	TH2F *fHistNPPPi0GammaEtaPtET;//Marcelo please add comment 
	TH2F *fHistNPPPi0GammaEtaET;//Marcelo please add comment 
	TH2F *fHistNPPPi0GammaEtaE;//Marcelo please add comment 
	TH2F *fHistNPPPi0GammaEtaPt;//Marcelo please add comment 
	TH1F *fHistNPPPi0GammatotET;//Marcelo please add comment 
		
	// *******************
	// electron ET inside EMCal acceptance
	// *******************
	TH2F *fHistElectronAccEtaEET;//Marcelo please add comment 
	TH2F *fHistElectronAccEtaPtET;//Marcelo please add comment 
	TH2F *fHistElectronAccEtaET;//Marcelo please add comment 
	TH2F *fHistElectronAccEtaE;//Marcelo please add comment 
	TH2F *fHistElectronAccEtaPt;//Marcelo please add comment 
	TH1F *fHistElectronAcctotET;//Marcelo please add comment 
	
	TH2F *fHistConvElectronAccEtaEET;//Marcelo please add comment  
	TH2F *fHistConvElectronAccEtaPtET;//Marcelo please add comment  
	TH2F *fHistConvElectronAccEtaET;//Marcelo please add comment  
	TH2F *fHistConvElectronAccEtaE;//Marcelo please add comment  
	TH2F *fHistConvElectronAccEtaPt;//Marcelo please add comment  
	TH1F *fHistConvElectronAcctotET;//Marcelo please add comment  
	
	TH2F *fHistScatElectronAccEtaEET;//Marcelo please add comment  
	TH2F *fHistScatElectronAccEtaPtET;//Marcelo please add comment  
	TH2F *fHistScatElectronAccEtaET;//Marcelo please add comment  
	TH2F *fHistScatElectronAccEtaE;//Marcelo please add comment  
	TH2F *fHistScatElectronAccEtaPt;//Marcelo please add comment  
	TH1F *fHistScatElectronAcctotET;//Marcelo please add comment  
	
	// *******************
	// total electron ET inside EMCal acceptance
	// *******************
	TH1F *fHistTotElectronAcctotET;//Marcelo please add comment

	// *******************
	// gamma ET inside EMCal acceptance
	// *******************
	TH2F *fHistGammaAccEtaEET;//Marcelo please add comment  
	TH2F *fHistGammaAccEtaPtET;//Marcelo please add comment  
	TH2F *fHistGammaAccEtaET;//Marcelo please add comment  
	TH2F *fHistGammaAccEtaE;//Marcelo please add comment  
	TH2F *fHistGammaAccEtaPt;//Marcelo please add comment  
	TH1F *fHistGammaAcctotET;//Marcelo please add comment  
	
	TH2F *fHistAnnihGammaAccEtaEET;//Marcelo please add comment  
	TH2F *fHistAnnihGammaAccEtaPtET;//Marcelo please add comment  
	TH2F *fHistAnnihGammaAccEtaET;//Marcelo please add comment  
	TH2F *fHistAnnihGammaAccEtaE;//Marcelo please add comment  
	TH2F *fHistAnnihGammaAccEtaPt;//Marcelo please add comment  
	TH1F *fHistAnnihGammaAcctotET;//Marcelo please add comment  
	
	TH2F *fHistScatGammaAccEtaEET;//Marcelo please add comment  
	TH2F *fHistScatGammaAccEtaPtET;//Marcelo please add comment  
	TH2F *fHistScatGammaAccEtaET;//Marcelo please add comment  
	TH2F *fHistScatGammaAccEtaE;//Marcelo please add comment  
	TH2F *fHistScatGammaAccEtaPt;//Marcelo please add comment  
	TH1F *fHistScatGammaAcctotET;//Marcelo please add comment  
	
	TH2F *fHistConvGammaAccEtaEET;//Marcelo please add comment  
	TH2F *fHistConvGammaAccEtaPtET;//Marcelo please add comment  
	TH2F *fHistConvGammaAccEtaET;//Marcelo please add comment  
	TH2F *fHistConvGammaAccEtaE;//Marcelo please add comment  
	TH2F *fHistConvGammaAccEtaPt;//Marcelo please add comment  
	TH1F *fHistConvGammaAcctotET;//Marcelo please add comment  
	
	TH2F *fHistNonConvGammaAccEtaEET;//Marcelo please add comment  
	TH2F *fHistNonConvGammaAccEtaPtET;//Marcelo please add comment  
	TH2F *fHistNonConvGammaAccEtaET;//Marcelo please add comment  
	TH2F *fHistNonConvGammaAccEtaE;//Marcelo please add comment  
	TH2F *fHistNonConvGammaAccEtaPt;//Marcelo please add comment  
	TH1F *fHistNonConvGammaAcctotET;//Marcelo please add comment  
	
	// *******************
	// total gamma ET inside EMCal acceptance
	// *******************
	TH1F *fHistTotGammaAcctotET;//Marcelo please add comment

	// *******************
	// total electromagnetic ET inside EMCal acceptance
	// *******************
	TH1F *fHistTotEMAcctotET;//Marcelo please add comment

	// non-primary electromagnetic ET
	TH2F *fHistNPPElectronAccEtaEET;//Marcelo please add comment 
	TH2F *fHistNPPElectronAccEtaPtET;//Marcelo please add comment 
	TH2F *fHistNPPElectronAccEtaE;//Marcelo please add comment 
	TH2F *fHistNPPElectronAccEtaPt;//Marcelo please add comment 
	
	TH2F *fHistNPPGammaAccEtaEET;//Marcelo please add comment 
	TH2F *fHistNPPGammaAccEtaPtET;//Marcelo please add comment 
	TH2F *fHistNPPGammaAccEtaE;//Marcelo please add comment 
	TH2F *fHistNPPGammaAccEtaPt;//Marcelo please add comment 	
	
	// *******************
	// electron ET reconstructed in EMCal
	// *******************
	TH2F *fHistElectronRecEtaEET;//Marcelo please add comment 
	TH2F *fHistElectronRecEtaPtET;//Marcelo please add comment 
	TH2F *fHistElectronRecEtaET;//Marcelo please add comment 
	TH2F *fHistElectronRecEtaE;//Marcelo please add comment 
	TH2F *fHistElectronRecEtaPt;//Marcelo please add comment 
	TH1F *fHistElectronRectotET;//Marcelo please add comment 
	
	TH2F *fHistConvElectronRecEtaEET;//Marcelo please add comment  
	TH2F *fHistConvElectronRecEtaPtET;//Marcelo please add comment  
	TH2F *fHistConvElectronRecEtaET;//Marcelo please add comment  
	TH2F *fHistConvElectronRecEtaE;//Marcelo please add comment  
	TH2F *fHistConvElectronRecEtaPt;//Marcelo please add comment  
	TH1F *fHistConvElectronRectotET;//Marcelo please add comment  
	
	TH2F *fHistScatElectronRecEtaEET;//Marcelo please add comment  
	TH2F *fHistScatElectronRecEtaPtET;//Marcelo please add comment  
	TH2F *fHistScatElectronRecEtaET;//Marcelo please add comment  
	TH2F *fHistScatElectronRecEtaE;//Marcelo please add comment  
	TH2F *fHistScatElectronRecEtaPt;//Marcelo please add comment  
	TH1F *fHistScatElectronRectotET;//Marcelo please add comment  
	
	// *******************
	// total Electron ET reconstructed in EMCal
	// *******************
	TH1F *fHistTotElectronRectotET;//Marcelo please add comment

	// *******************
	// gamma ET reconstructed in EMCal
	// *******************
	TH2F *fHistGammaRecEtaEET;//Marcelo please add comment  
	TH2F *fHistGammaRecEtaPtET;//Marcelo please add comment  
	TH2F *fHistGammaRecEtaET;//Marcelo please add comment  
	TH2F *fHistGammaRecEtaE;//Marcelo please add comment  
	TH2F *fHistGammaRecEtaPt;//Marcelo please add comment  
	TH1F *fHistGammaRectotET;//Marcelo please add comment  
	
	TH2F *fHistAnnihGammaRecEtaEET;//Marcelo please add comment  
	TH2F *fHistAnnihGammaRecEtaPtET;//Marcelo please add comment  
	TH2F *fHistAnnihGammaRecEtaET;//Marcelo please add comment  
	TH2F *fHistAnnihGammaRecEtaE;//Marcelo please add comment  
	TH2F *fHistAnnihGammaRecEtaPt;//Marcelo please add comment  
	TH1F *fHistAnnihGammaRectotET;//Marcelo please add comment  
	
	TH2F *fHistScatGammaRecEtaEET;//Marcelo please add comment  
	TH2F *fHistScatGammaRecEtaPtET;//Marcelo please add comment  
	TH2F *fHistScatGammaRecEtaET;//Marcelo please add comment  
	TH2F *fHistScatGammaRecEtaE;//Marcelo please add comment  
	TH2F *fHistScatGammaRecEtaPt;//Marcelo please add comment  
	TH1F *fHistScatGammaRectotET;//Marcelo please add comment  

	// *******************
	// total gamma ET reconstructed in EMCal
	// *******************
	TH1F *fHistTotGammaRectotET;//Marcelo please add comment

	// *******************
	// total EM ET reconstructed in EMCal
	// *******************
	TH1F *fHistTotEMRectotET;//Marcelo please add comment

	// non-primary electromagnetic ET
	TH2F *fHistNPPElectronRecEtaEET;//Marcelo please add comment 
	TH2F *fHistNPPElectronRecEtaPtET;//Marcelo please add comment 
	TH2F *fHistNPPElectronRecEtaET;//Marcelo please add comment 
	TH2F *fHistNPPElectronRecEtaE;//Marcelo please add comment 
	TH2F *fHistNPPElectronRecEtaPt;//Marcelo please add comment 
	TH1F *fHistNPPElectronRectotET;//Marcelo please add comment 
	
	TH2F *fHistNPPGammaRecEtaEET;//Marcelo please add comment 
	TH2F *fHistNPPGammaRecEtaPtET;//Marcelo please add comment 
	TH2F *fHistNPPGammaRecEtaET;//Marcelo please add comment 
	TH2F *fHistNPPGammaRecEtaE;//Marcelo please add comment 
	TH2F *fHistNPPGammaRecEtaPt;//Marcelo please add comment 
	TH1F *fHistNPPGammaRectotET;//Marcelo please add comment 
	
	TH1F *fHistTotNPPEMRectotET;//Marcelo please add comment

	TH2F *fHistNPPPi0GammaRecEtaEET;//Marcelo please add comment 
	TH2F *fHistNPPPi0GammaRecEtaPtET;//Marcelo please add comment 
	TH2F *fHistNPPPi0GammaRecEtaET;//Marcelo please add comment 
	TH2F *fHistNPPPi0GammaRecEtaE;//Marcelo please add comment 
	TH2F *fHistNPPPi0GammaRecEtaPt;//Marcelo please add comment 
	TH1F *fHistNPPPi0GammaRectotET;//Marcelo please add comment 
	
	// *******************
	// muon ET (+ and -)
	// *******************
	TH2F *fHistMuonEtaEET;//Marcelo please add comment 
	TH2F *fHistMuonAccEtaEET;//Marcelo please add comment 
	TH2F *fHistMuonRecEtaEET;//Marcelo please add comment 
	TH2F *fHistMuonMatchEtaEET;//Marcelo please add comment 

	TH2F *fHistMuonEtaPtET;//Marcelo please add comment 
	TH2F *fHistMuonAccEtaPtET;//Marcelo please add comment 
	TH2F *fHistMuonRecEtaPtET;//Marcelo please add comment 
	TH2F *fHistMuonMatchEtaPtET;//Marcelo please add comment 

	TH2F *fHistMuonEtaET;//Marcelo please add comment 
	TH2F *fHistMuonAccEtaET;//Marcelo please add comment 
	TH2F *fHistMuonRecEtaET;//Marcelo please add comment 
	TH2F *fHistMuonMatchEtaET;//Marcelo please add comment 
	
	TH2F *fHistMuonEtaE;//Marcelo please add comment 
	TH2F *fHistMuonAccEtaE;//Marcelo please add comment 
	TH2F *fHistMuonRecEtaE;//Marcelo please add comment 
	TH2F *fHistMuonMatchEtaE;//Marcelo please add comment 
	
	TH2F *fHistMuonEtaPt;//Marcelo please add comment 
	TH2F *fHistMuonAccEtaPt;//Marcelo please add comment 
	TH2F *fHistMuonRecEtaPt;//Marcelo please add comment 
	TH2F *fHistMuonMatchEtaPt;//Marcelo please add comment 
	
	TH1F *fHistMuontotET;//Marcelo please add comment 
	TH1F *fHistMuonAcctotET;//Marcelo please add comment 
	TH1F *fHistMuonRectotET;//Marcelo please add comment 
	TH1F *fHistMuonMatchtotET;//Marcelo please add comment 
	
	TH1F *fHistMuonRectotETDep;//Marcelo please add comment 
	TH1F *fHistMuonMatchtotETDep;//Marcelo please add comment 
	
	TH2F *fHistMuonRecEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistMuonMatchEtaEDepETDep;//Marcelo please add comment 

	TH2F *fHistMuonRecEtaPtETDep;//Marcelo please add comment 
	TH2F *fHistMuonMatchEtaPtETDep;//Marcelo please add comment 
	
	TH2F *fHistMuonRecEtaETDep;//Marcelo please add comment 
	TH2F *fHistMuonMatchEtaETDep;//Marcelo please add comment 

	TH2F *fHistMuonRecResEET;//Marcelo please add comment 
	TH2F *fHistMuonRecResPtET;//Marcelo please add comment 
	TH2F *fHistMuonRecResE;//Marcelo please add comment 
	TH2F *fHistMuonRecResPt;//Marcelo please add comment 
	TH2F *fHistMuonRecResEDepETDep;//Marcelo please add comment 
	TH2F *fHistMuonRecResPtETDep;//Marcelo please add comment 
	
	// *******************
	// pion ET (+ and -)
	// *******************
	TH2F *fHistPionEtaEET;//Marcelo please add comment 
	TH2F *fHistPionAccEtaEET;//Marcelo please add comment 
	TH2F *fHistPionRecEtaEET;//Marcelo please add comment 
	TH2F *fHistPionMatchEtaEET;//Marcelo please add comment 
	
	TH2F *fHistPionEtaPtET;//Marcelo please add comment 
	TH2F *fHistPionAccEtaPtET;//Marcelo please add comment 
	TH2F *fHistPionRecEtaPtET;//Marcelo please add comment 
	TH2F *fHistPionMatchEtaPtET;//Marcelo please add comment 
	
	TH2F *fHistPionEtaET;//Marcelo please add comment 
	TH2F *fHistPionAccEtaET;//Marcelo please add comment 
	TH2F *fHistPionRecEtaET;//Marcelo please add comment 
	TH2F *fHistPionMatchEtaET;//Marcelo please add comment 
	
	TH2F *fHistPionEtaE;//Marcelo please add comment 
	TH2F *fHistPionAccEtaE;//Marcelo please add comment 
	TH2F *fHistPionRecEtaE;//Marcelo please add comment 
	TH2F *fHistPionMatchEtaE;//Marcelo please add comment 
	
	TH2F *fHistPionEtaPt;//Marcelo please add comment 
	TH2F *fHistPionAccEtaPt;//Marcelo please add comment 
	TH2F *fHistPionRecEtaPt;//Marcelo please add comment 
	TH2F *fHistPionMatchEtaPt;//Marcelo please add comment 
	
	TH1F *fHistPiontotET;//Marcelo please add comment 
	TH1F *fHistPionAcctotET;//Marcelo please add comment 
	TH1F *fHistPionRectotET;//Marcelo please add comment 
	TH1F *fHistPionMatchtotET;//Marcelo please add comment 
	
	TH1F *fHistPionRectotETDep;//Marcelo please add comment 
	TH1F *fHistPionMatchtotETDep;//Marcelo please add comment 
	
	TH2F *fHistPionRecEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistPionMatchEtaEDepETDep;//Marcelo please add comment 

	TH2F *fHistPionRecEtaPtETDep;//Marcelo please add comment 
	TH2F *fHistPionMatchEtaPtETDep;//Marcelo please add comment 
	
	TH2F *fHistPionRecEtaETDep;//Marcelo please add comment 
	TH2F *fHistPionMatchEtaETDep;//Marcelo please add comment 
	
	TH2F *fHistPionRecResEET;//Marcelo please add comment 
	TH2F *fHistPionRecResPtET;//Marcelo please add comment 
	TH2F *fHistPionRecResE;//Marcelo please add comment 
	TH2F *fHistPionRecResPt;//Marcelo please add comment 
	TH2F *fHistPionRecResEDepETDep;//Marcelo please add comment 
	TH2F *fHistPionRecResPtETDep;//Marcelo please add comment 
	
	// *******************
	// charged kaon (+ and -) ET
	// *******************
	TH2F *fHistKaonEtaEET;//Marcelo please add comment 
	TH2F *fHistKaonAccEtaEET;//Marcelo please add comment 
	TH2F *fHistKaonRecEtaEET;//Marcelo please add comment 
	TH2F *fHistKaonMatchEtaEET;//Marcelo please add comment 
	
	TH2F *fHistKaonEtaPtET;//Marcelo please add comment 
	TH2F *fHistKaonAccEtaPtET;//Marcelo please add comment 
	TH2F *fHistKaonRecEtaPtET;//Marcelo please add comment 
	TH2F *fHistKaonMatchEtaPtET;//Marcelo please add comment 
	
	TH2F *fHistKaonEtaET;//Marcelo please add comment 
	TH2F *fHistKaonAccEtaET;//Marcelo please add comment 
	TH2F *fHistKaonRecEtaET;//Marcelo please add comment 
	TH2F *fHistKaonMatchEtaET;//Marcelo please add comment 
	
	TH2F *fHistKaonEtaE;//Marcelo please add comment 
	TH2F *fHistKaonAccEtaE;//Marcelo please add comment 
	TH2F *fHistKaonRecEtaE;//Marcelo please add comment 
	TH2F *fHistKaonMatchEtaE;//Marcelo please add comment 
	
	TH2F *fHistKaonEtaPt;//Marcelo please add comment 
	TH2F *fHistKaonAccEtaPt;//Marcelo please add comment 
	TH2F *fHistKaonRecEtaPt;//Marcelo please add comment 
	TH2F *fHistKaonMatchEtaPt;//Marcelo please add comment 

	TH1F *fHistKaontotET;//Marcelo please add comment 
	TH1F *fHistKaonAcctotET;//Marcelo please add comment 
	TH1F *fHistKaonRectotET;//Marcelo please add comment 
	TH1F *fHistKaonMatchtotET;//Marcelo please add comment 
	
	TH1F *fHistKaonRectotETDep;//Marcelo please add comment 
	TH1F *fHistKaonMatchtotETDep;//Marcelo please add comment 
	
	TH2F *fHistKaonRecEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistKaonMatchEtaEDepETDep;//Marcelo please add comment 

	TH2F *fHistKaonRecEtaPtETDep;//Marcelo please add comment 
	TH2F *fHistKaonMatchEtaPtETDep;//Marcelo please add comment 
	
	TH2F *fHistKaonRecEtaETDep;//Marcelo please add comment 
	TH2F *fHistKaonMatchEtaETDep;//Marcelo please add comment 
	
	TH2F *fHistKaonRecResEET;//Marcelo please add comment 
	TH2F *fHistKaonRecResPtET;//Marcelo please add comment 
	TH2F *fHistKaonRecResE;//Marcelo please add comment 
	TH2F *fHistKaonRecResPt;//Marcelo please add comment 
	TH2F *fHistKaonRecResEDepETDep;//Marcelo please add comment 
	TH2F *fHistKaonRecResPtETDep;//Marcelo please add comment 	
	
	// *******************
	// proton (anti) ET
	// *******************
	TH2F *fHistProtonEtaEET;//Marcelo please add comment 
	TH2F *fHistProtonAccEtaEET;//Marcelo please add comment 
	TH2F *fHistProtonRecEtaEET;//Marcelo please add comment 
	TH2F *fHistProtonMatchEtaEET;//Marcelo please add comment 
	
	TH2F *fHistProtonEtaPtET;//Marcelo please add comment 
	TH2F *fHistProtonAccEtaPtET;//Marcelo please add comment 
	TH2F *fHistProtonRecEtaPtET;//Marcelo please add comment 
	TH2F *fHistProtonMatchEtaPtET;//Marcelo please add comment 
	
	TH2F *fHistProtonEtaET;//Marcelo please add comment 
	TH2F *fHistProtonAccEtaET;//Marcelo please add comment 
	TH2F *fHistProtonRecEtaET;//Marcelo please add comment 
	TH2F *fHistProtonMatchEtaET;//Marcelo please add comment 
	
	TH2F *fHistProtonEtaE;//Marcelo please add comment 
	TH2F *fHistProtonAccEtaE;//Marcelo please add comment 
	TH2F *fHistProtonRecEtaE;//Marcelo please add comment 
	TH2F *fHistProtonMatchEtaE;//Marcelo please add comment 
	
	TH2F *fHistProtonEtaPt;//Marcelo please add comment 
	TH2F *fHistProtonAccEtaPt;//Marcelo please add comment 
	TH2F *fHistProtonRecEtaPt;//Marcelo please add comment 
	TH2F *fHistProtonMatchEtaPt;//Marcelo please add comment 

	TH1F *fHistProtontotET;//Marcelo please add comment 
	TH1F *fHistProtonAcctotET;//Marcelo please add comment 
	TH1F *fHistProtonRectotET;//Marcelo please add comment 
	TH1F *fHistProtonMatchtotET;//Marcelo please add comment 
	
	TH1F *fHistProtonRectotETDep;//Marcelo please add comment 
	TH1F *fHistProtonMatchtotETDep;//Marcelo please add comment 
	
	TH2F *fHistProtonRecEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistProtonMatchEtaEDepETDep;//Marcelo please add comment 
	
	TH2F *fHistProtonRecEtaPtETDep;//Marcelo please add comment 
	TH2F *fHistProtonMatchEtaPtETDep;//Marcelo please add comment 
	
	TH2F *fHistProtonRecEtaETDep;//Marcelo please add comment 
	TH2F *fHistProtonMatchEtaETDep;//Marcelo please add comment 

	TH2F *fHistProtonRecResEET;//Marcelo please add comment 
	TH2F *fHistProtonRecResPtET;//Marcelo please add comment 
	TH2F *fHistProtonRecResE;//Marcelo please add comment 
	TH2F *fHistProtonRecResPt;//Marcelo please add comment 
	TH2F *fHistProtonRecResEDepETDep;//Marcelo please add comment 
	TH2F *fHistProtonRecResPtETDep;//Marcelo please add comment 
	
	// *******************
	// total charged ET
	// *******************
	TH1F *fHistTotChargedtotET;//Marcelo please add comment
	TH1F *fHistTotChargedAcctotET;//Marcelo please add comment
	TH1F *fHistTotChargedRectotET;//Marcelo please add comment
	TH1F *fHistTotChargedRectotETDep;//Marcelo please add comment
	TH1F *fHistTotChargedMatchtotET;//Marcelo please add comment
	TH1F *fHistTotChargedMatchtotETDep;//Marcelo please add comment
	
	// *******************
	// neutron (anti) ET
	// *******************
	TH2F *fHistNeutronEtaEET;//Marcelo please add comment 
	TH2F *fHistNeutronAccEtaEET;//Marcelo please add comment 
	TH2F *fHistNeutronRecEtaEET;//Marcelo please add comment 
	
	TH2F *fHistNeutronEtaPtET;//Marcelo please add comment 
	TH2F *fHistNeutronAccEtaPtET;//Marcelo please add comment 
	TH2F *fHistNeutronRecEtaPtET;//Marcelo please add comment 
	
	TH2F *fHistNeutronEtaET;//Marcelo please add comment 
	TH2F *fHistNeutronAccEtaET;//Marcelo please add comment 
	TH2F *fHistNeutronRecEtaET;//Marcelo please add comment 
	
	TH2F *fHistNeutronEtaE;//Marcelo please add comment 
	TH2F *fHistNeutronAccEtaE;//Marcelo please add comment 
	TH2F *fHistNeutronRecEtaE;//Marcelo please add comment 
	
	TH2F *fHistNeutronEtaPt;//Marcelo please add comment 
	TH2F *fHistNeutronAccEtaPt;//Marcelo please add comment 
	TH2F *fHistNeutronRecEtaPt;//Marcelo please add comment 
	
	TH1F *fHistNeutrontotET;//Marcelo please add comment 
	TH1F *fHistNeutronAcctotET;//Marcelo please add comment 
	TH1F *fHistNeutronRectotET;//Marcelo please add comment 
	TH1F *fHistNeutronRectotETDep;//Marcelo please add comment 
	
	TH2F *fHistNeutronRecEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistNeutronRecEtaETDep;//Marcelo please add comment 
	
	TH2F *fHistNeutronRecEtaPtETDep;//Marcelo please add comment 
		
	// *******************
	// neutral kaon ET
	// *******************
	TH2F *fHistK0EtaEET;//Marcelo please add comment 
	TH2F *fHistK0RecEtaEET;//Marcelo please add comment 
	
	TH2F *fHistK0EtaPtET;//Marcelo please add comment 
	TH2F *fHistK0RecEtaPtET;//Marcelo please add comment 
	
	TH2F *fHistK0EtaET;//Marcelo please add comment 
	TH2F *fHistK0RecEtaET;//Marcelo please add comment 
	
	TH2F *fHistK0EtaE;//Marcelo please add comment 
	TH2F *fHistK0RecEtaE;//Marcelo please add comment 
	
	TH2F *fHistK0EtaPt;//Marcelo please add comment 
	TH2F *fHistK0RecEtaPt;//Marcelo please add comment 

	TH1F *fHistK0totET;//Marcelo please add comment 
	TH1F *fHistK0RectotET;//Marcelo please add comment 
	
	TH1F *fHistK0RectotETDep;//Marcelo please add comment 
	
	TH2F *fHistK0RecEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistK0RecEtaETDep;//Marcelo please add comment 
	
	TH2F *fHistK0RecEtaPtETDep;//Marcelo please add comment 
		
	// *******************
	// Lambda(anti) ET
	// *******************
	TH2F *fHistLambdaEtaEET;//Marcelo please add comment 
	TH2F *fHistLambdaRecEtaEET;//Marcelo please add comment 
	
	TH2F *fHistLambdaEtaPtET;//Marcelo please add comment 
	TH2F *fHistLambdaRecEtaPtET;//Marcelo please add comment 
	
	TH2F *fHistLambdaEtaET;//Marcelo please add comment 
	TH2F *fHistLambdaRecEtaET;//Marcelo please add comment 
	
	TH2F *fHistLambdaEtaE;//Marcelo please add comment 
	TH2F *fHistLambdaRecEtaE;//Marcelo please add comment 
	
	TH2F *fHistLambdaEtaPt;//Marcelo please add comment 
	TH2F *fHistLambdaRecEtaPt;//Marcelo please add comment 
	
	TH1F *fHistLambdatotET;//Marcelo please add comment 
	TH1F *fHistLambdaRectotET;//Marcelo please add comment 
	
	TH1F *fHistLambdaRectotETDep;//Marcelo please add comment 
	
	TH2F *fHistLambdaRecEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistLambdaRecEtaETDep;//Marcelo please add comment 
	
	TH2F *fHistLambdaRecEtaPtETDep;//Marcelo please add comment 

	// *******************
	// total neutral ET
	// *******************
	TH1F *fHistTotNeutraltotET;//Marcelo please add comment
	TH1F *fHistTotNeutralRectotET;//Marcelo please add comment
	TH1F *fHistTotNeutralRectotETDep;//Marcelo please add comment
	
	// *******************
	// total ET
	// *******************
	TH1F *fHistTotaltotET;//Marcelo please add comment
	TH1F *fHistTotalAcctotET;//Marcelo please add comment
	TH1F *fHistTotalRectotET;//Marcelo please add comment
	TH1F *fHistTotalRectotETDep;//Marcelo please add comment
	
	// *******************
	// some checks
	// *******************

	// check produced electrons
	//TH1F *fHistElectronFirstMother;//Marcelo please add comment 
	TH2F *fHistElectronFirstMotherXY;//Marcelo please add comment 
	TH1F *fHistElectronNDaughters;//Marcelo please add comment 
	TH1F *fHistElectronDaughters;//Marcelo please add comment 
	TH2F *fHistElectronDaughtersXY;//Marcelo please add comment 

	TH1F *fHistElectronFirstMotherAcc;//Marcelo please add comment  
	TH2F *fHistElectronFirstMotherXYAcc;//Marcelo please add comment  
	TH1F *fHistElectronNDaughtersAcc;//Marcelo please add comment 
	TH1F *fHistElectronDaughtersAcc;//Marcelo please add comment 
	TH2F *fHistElectronDaughtersXYAcc;//Marcelo please add comment 

	TH1F *fHistElectronFirstMotherRec;//Marcelo please add comment  
	TH2F *fHistElectronFirstMotherXYRec;//Marcelo please add comment  
	TH1F *fHistElectronNDaughtersRec;//Marcelo please add comment 
	TH1F *fHistElectronDaughtersRec;//Marcelo please add comment 
	TH2F *fHistElectronDaughtersXYRec;//Marcelo please add comment 

	TH1F *fHistNPPElectronFirstMother;//Marcelo please add comment 
	TH2F *fHistNPPElectronFirstMotherXY;//Marcelo please add comment 
	TH1F *fHistNPPElectronNDaughters;//Marcelo please add comment 
	TH1F *fHistNPPElectronDaughters;//Marcelo please add comment 
	TH2F *fHistNPPElectronDaughtersXY;//Marcelo please add comment 
	
	TH1F *fHistNPPElectronFirstMotherAcc;//Marcelo please add comment  
	TH2F *fHistNPPElectronFirstMotherXYAcc;//Marcelo please add comment  
	TH1F *fHistNPPElectronNDaughtersAcc;//Marcelo please add comment 
	TH1F *fHistNPPElectronDaughtersAcc;//Marcelo please add comment 
	TH2F *fHistNPPElectronDaughtersXYAcc;//Marcelo please add comment 
	
	TH1F *fHistNPPElectronFirstMotherRec;//Marcelo please add comment  
	TH2F *fHistNPPElectronFirstMotherXYRec;//Marcelo please add comment  
	TH1F *fHistNPPElectronNDaughtersRec;//Marcelo please add comment 
	TH1F *fHistNPPElectronDaughtersRec;//Marcelo please add comment 
	TH2F *fHistNPPElectronDaughtersXYRec;//Marcelo please add comment 
	
	// check produced gammas
	//TH1F *fHistGammaFirstMother;//Marcelo please add comment 
	TH2F *fHistGammaFirstMotherXY;//Marcelo please add comment 
	TH1F *fHistGammaNDaughters;//Marcelo please add comment 
	TH1F *fHistGammaDaughters;//Marcelo please add comment 
	TH2F *fHistGammaDaughtersXY;//Marcelo please add comment 
	TH2F *fHistConvGammaDaughtersXY;//Marcelo please add comment 
	TH2F *fHistNonConvGammaDaughtersXY;//Marcelo please add comment 
	
	TH1F *fHistGammaFirstMotherAcc;//Marcelo please add comment  
	TH2F *fHistGammaFirstMotherXYAcc;//Marcelo please add comment  
	TH1F *fHistGammaNDaughtersAcc;//Marcelo please add comment 
	TH1F *fHistGammaDaughtersAcc;//Marcelo please add comment 
	TH2F *fHistGammaDaughtersXYAcc;//Marcelo please add comment 
	TH2F *fHistConvGammaDaughtersXYAcc;//Marcelo please add comment 
	TH2F *fHistNonConvGammaDaughtersXYAcc;//Marcelo please add comment 
	
	TH1F *fHistGammaFirstMotherRec;//Marcelo please add comment  
	TH2F *fHistGammaFirstMotherXYRec;//Marcelo please add comment  
	TH1F *fHistGammaNDaughtersRec;//Marcelo please add comment 
	TH1F *fHistGammaDaughtersRec;//Marcelo please add comment 
	TH2F *fHistGammaDaughtersXYRec;//Marcelo please add comment 
	TH2F *fHistConvGammaDaughtersXYRec;//Marcelo please add comment 
	TH2F *fHistNonConvGammaDaughtersXYRec;//Marcelo please add comment 
	
	TH1F *fHistNPPGammaFirstMother;//Marcelo please add comment 
	TH2F *fHistNPPGammaFirstMotherXY;//Marcelo please add comment 
	TH1F *fHistNPPGammaNDaughters;//Marcelo please add comment 
	TH1F *fHistNPPGammaDaughters;//Marcelo please add comment 
	TH2F *fHistNPPGammaDaughtersXY;//Marcelo please add comment 
	
	TH1F *fHistNPPGammaFirstMotherAcc;//Marcelo please add comment  
	TH2F *fHistNPPGammaFirstMotherXYAcc;//Marcelo please add comment  
	TH1F *fHistNPPGammaNDaughtersAcc;//Marcelo please add comment 
	TH1F *fHistNPPGammaDaughtersAcc;//Marcelo please add comment 
	TH2F *fHistNPPGammaDaughtersXYAcc;//Marcelo please add comment 
	
	TH1F *fHistNPPGammaFirstMotherRec;//Marcelo please add comment  
	TH2F *fHistNPPGammaFirstMotherXYRec;//Marcelo please add comment  
	TH1F *fHistNPPGammaNDaughtersRec;//Marcelo please add comment 
	TH1F *fHistNPPGammaDaughtersRec;//Marcelo please add comment 
	TH2F *fHistNPPGammaDaughtersXYRec;//Marcelo please add comment 

	//check projections
	//TH2F *fHistAllERecEMC;//Marcelo please add comment	
	TH2F *fHistAllPtRecPtMC;//Marcelo please add comment
	//TH2F *fHistElectronERecEMC;//Marcelo please add comment	
	//TH2F *fHistGammaERecEMC;//Marcelo please add comment
	
	TH2F *fHistChargedRes;//Marcelo please add comment
	TH2F *fHistChargedRes2;//Marcelo please add comment
	TH2F *fHistChargedRes3;//Marcelo please add comment
	TH2F *fHistNeutralRes;//Marcelo please add comment
	TH2F *fHistElectronRes;//Marcelo please add comment
	TH2F *fHistGammaRes;//Marcelo please add comment
	
	TH2F *fHistIsInAcc;//Marcelo please add comment
	
 private:

  //Declare it private to avoid compilation warning
    AliAnalysisEmEtMonteCarlo & operator = (const AliAnalysisEmEtMonteCarlo & g) ;//cpy assignment
    AliAnalysisEmEtMonteCarlo(const AliAnalysisEmEtMonteCarlo & g) ; // cpy ctor
    ClassDef(AliAnalysisEmEtMonteCarlo, 1);
};

#endif //ALIANALYSISEMETMONTECARLO_H
