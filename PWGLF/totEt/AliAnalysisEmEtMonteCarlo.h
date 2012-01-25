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

	Double_t fResCut; // b(fm), for Hijing; 0 otherwise
	Double_t fPrimtotET, fPrimAcctotET, fPrimRectotET, fPrimRectotETDep;//Primary particles ET

	Double_t fElectrontotET, fElectronAcctotET, fElectronRectotET;//Electron ET
	Double_t fConvElectrontotET, fConvElectronAcctotET, fConvElectronRectotET, fScatElectrontotET, fScatElectronAcctotET, fScatElectronRectotET;//Secondary electrons ET
	Double_t fTotElectrontotET, fTotElectronAcctotET, fTotElectronRectotET;//Total Electrons ET

	Double_t fGammatotET, fGammaAcctotET, fGammaRectotET;//Gamma event ET
	Double_t fAnnihGammatotET, fAnnihGammaAcctotET, fAnnihGammaRectotET, fScatGammatotET, fScatGammaAcctotET, fScatGammaRectotET;//Secondary gamma ET
	Double_t fTotGammatotET, fTotGammaAcctotET, fTotGammaRectotET;//Total gamma ET
	Double_t fConvGammatotET, fNonConvGammatotET, fConvGammaAcctotET, fNonConvGammaAcctotET, fNPPPi0GammatotET, fNPPPi0GammaRectotET;//Conversion and pi0 gamma ET

	Double_t fTotEMtotET, fTotEMAcctotET, fTotEMRectotET;//Total EM ET

	Double_t fNPPElectrontotET, fNPPElectronRectotET, fNPPGammatotET, fNPPGammaRectotET;//Non-physical primary electron ET
	Double_t fTotNPPEMtotET, fTotNPPEMRectotET;//Total Non-physical primary electron ET

	Double_t fMuontotET, fPiontotET, fKaontotET, fProtontotET;//Charged particles ET
	Double_t fMuonAcctotET, fPionAcctotET, fKaonAcctotET, fProtonAcctotET;//Charged particles acceptance ET
	Double_t fMuonRectotET, fMuonRectotETDep, fPionRectotET, fPionRectotETDep, fKaonRectotET, fKaonRectotETDep, fProtonRectotET, fProtonRectotETDep;//Charged particles reconstructed ET
	Double_t fMuonMatchtotET, fMuonMatchtotETDep, fPionMatchtotET, fPionMatchtotETDep, fKaonMatchtotET, fKaonMatchtotETDep, fProtonMatchtotET, fProtonMatchtotETDep;//Charged particles track matched ET
	Double_t fTotChargedtotET, fTotChargedAcctotET, fTotChargedRectotET, fTotChargedRectotETDep, fTotChargedMatchtotET, fTotChargedMatchtotETDep;//Total charged particles ET

	Double_t fNeutrontotET, fNeutronAcctotET, fNeutronRectotET, fNeutronRectotETDep;//Neutrons ET
	Double_t fK0totET, fK0RectotET, fK0RectotETDep, fLambdatotET, fLambdaRectotET, fLambdaRectotETDep;//K0 and Lambda ET
	Double_t fTotNeutraltotET, fTotNeutralRectotET, fTotNeutralRectotETDep;//Neutral particles ET

	Double_t fTotaltotET, fTotalAcctotET, fTotalRectotET, fTotalRectotETDep;//Total ET
	
	AliEMCALGeometry *fGeoUt;//EMCal geometry object

	// *******************
	// primaries ET
	// *******************
	TH2F *fHistPrimEtaEET;//total ET - Eta vs E 
	TH2F *fHistPrimEtaPtET;//total ET - Eta vs pt 
	TH2F *fHistPrimEtaET;//total ET - Eta
	TH1F *fHistPrimtotET;//total ET distribution
	
	TH2F *fHistPrimAccEtaEET;//acceptance ET - Eta vs E 
	TH2F *fHistPrimAccEtaPtET;//acceptance ET - Eta vs pt 
	TH2F *fHistPrimAccEtaET;//acceptance ET - Eta
	TH1F *fHistPrimAcctotET;//acceptance ET distribution
	
	TH2F *fHistPrimRecEtaEET;//reconstructed ET - Eta vs E 
	TH2F *fHistPrimRecEtaPtET;//reconstructed ET - Eta vs pt 
	TH2F *fHistPrimRecEtaET;//reconstructed ET - Eta
	TH1F *fHistPrimRectotET;//reconstructed ET distribution

	TH2F *fHistPrimRecEtaEDepETDep;//deposited ET - Eta vs E deposited
	TH2F *fHistPrimRecEtaPtETDep;//deposited ET - Eta vs pt 
	TH2F *fHistPrimRecEtaETDep;//deposited ET - Eta 
	TH1F *fHistPrimRectotETDep;//deposited ET distribution
	
	// *******************
	// electron ET
	// *******************
	TH2F *fHistElectronEtaEET;// ET - Eta vs E 
	TH2F *fHistElectronEtaPtET;//ET - Eta vs pt
	TH2F *fHistElectronEtaET;// ET - Eta 
	TH2F *fHistElectronEtaE;// multiplicity - Eta vs E
	TH2F *fHistElectronEtaPt;// multiplicity - Eta vs pt
	TH1F *fHistElectrontotET;// total ET distribution 

	TH2F *fHistConvElectronEtaEET;//ET - Eta vs E   
	TH2F *fHistConvElectronEtaPtET;//ET - Eta vs pt 
	TH2F *fHistConvElectronEtaET;//  ET - Eta
	TH2F *fHistConvElectronEtaE;//  multiplicity - Eta vs E
	TH2F *fHistConvElectronEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistConvElectrontotET;// total ET distribution 

	TH2F *fHistScatElectronEtaEET;//ET - Eta vs E  
	TH2F *fHistScatElectronEtaPtET;//ET - Eta vs pt 
	TH2F *fHistScatElectronEtaET;//  ET - Eta
	TH2F *fHistScatElectronEtaE;//  multiplicity - Eta vs E
	TH2F *fHistScatElectronEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistScatElectrontotET;//  total ET distribution
	
	// *******************
	// total electron ET
	// *******************
	TH1F *fHistTotElectrontotET;//total ET distribution
	
	// *******************
	// gamma ET
	// *******************
	TH2F *fHistGammaEtaEET;//ET - Eta vs E  
	TH2F *fHistGammaEtaPtET;//ET - Eta vs pt  
	TH2F *fHistGammaEtaET;//  ET - Eta
	TH2F *fHistGammaEtaE;//  multiplicity - Eta vs E
	TH2F *fHistGammaEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistGammatotET;//  total ET distribution
	
	TH2F *fHistAnnihGammaEtaEET;//ET - Eta vs E  
	TH2F *fHistAnnihGammaEtaPtET;//ET - Eta vs pt
	TH2F *fHistAnnihGammaEtaET;//  ET - Eta
	TH2F *fHistAnnihGammaEtaE;//  multiplicity - Eta vs E
	TH2F *fHistAnnihGammaEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistAnnihGammatotET;//  total ET distribution

	TH2F *fHistScatGammaEtaEET;//ET - Eta vs E   
	TH2F *fHistScatGammaEtaPtET;//ET - Eta vs pt  
	TH2F *fHistScatGammaEtaET;//  ET - Eta
	TH2F *fHistScatGammaEtaE;//  multiplicity - Eta vs E
	TH2F *fHistScatGammaEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistScatGammatotET;//  total ET distribution

	TH2F *fHistConvGammaEtaEET;//ET - Eta vs E  
	TH2F *fHistConvGammaEtaPtET;//ET - Eta vs pt
	TH2F *fHistConvGammaEtaET;//  ET - Eta
	TH2F *fHistConvGammaEtaE;//  multiplicity - Eta vs E
	TH2F *fHistConvGammaEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistConvGammatotET;//  total ET distribution
	
	TH2F *fHistNonConvGammaEtaEET;//ET - Eta vs E  
	TH2F *fHistNonConvGammaEtaPtET;//ET - Eta vs pt 
	TH2F *fHistNonConvGammaEtaET;//  ET - Eta
	TH2F *fHistNonConvGammaEtaE;//  multiplicity - Eta vs E
	TH2F *fHistNonConvGammaEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistNonConvGammatotET;//  total ET distribution
	
	// *******************
	// total gamma ET
	// *******************
	TH1F *fHistTotGammatotET;//total ET distribution

	// *******************
	// total electromagnetic ET
	// *******************
	TH1F *fHistTotEMtotET;//total ET distribution

	// non-primary electromagnetic ET
	TH2F *fHistNPPElectronEtaEET;//ET - Eta vs E 
	TH2F *fHistNPPElectronEtaPtET;//ET - Eta vs pt
	TH2F *fHistNPPElectronEtaET;// ET - Eta
	TH2F *fHistNPPElectronEtaE;// multiplicity - Eta vs E
	TH2F *fHistNPPElectronEtaPt;// multiplicity - Eta vs pt
	TH1F *fHistNPPElectrontotET;// total ET distribution

	TH2F *fHistNPPGammaEtaEET;//ET - Eta vs E  
	TH2F *fHistNPPGammaEtaPtET;//ET - Eta vs pt 
	TH2F *fHistNPPGammaEtaET;// ET - Eta
	TH2F *fHistNPPGammaEtaE;// multiplicity - Eta vs E
	TH2F *fHistNPPGammaEtaPt;// multiplicity - Eta vs pt
	TH1F *fHistNPPGammatotET;// total ET distribution

	TH1F *fHistTotNPPEMtotET;//total ET distribution

	TH2F *fHistNPPPi0GammaEtaEET;//ET - Eta vs E  
	TH2F *fHistNPPPi0GammaEtaPtET;//ET - Eta vs pt
	TH2F *fHistNPPPi0GammaEtaET;// ET - Eta
	TH2F *fHistNPPPi0GammaEtaE;// multiplicity - Eta vs E
	TH2F *fHistNPPPi0GammaEtaPt;// multiplicity - Eta vs pt
	TH1F *fHistNPPPi0GammatotET;// total ET distribution
		
	// *******************
	// electron ET inside EMCal acceptance
	// *******************
	TH2F *fHistElectronAccEtaEET;//ET - Eta vs E  
	TH2F *fHistElectronAccEtaPtET;//ET - Eta vs pt 
	TH2F *fHistElectronAccEtaET;// ET - Eta
	TH2F *fHistElectronAccEtaE;// multiplicity - Eta vs E
	TH2F *fHistElectronAccEtaPt;// multiplicity - Eta vs pt
	TH1F *fHistElectronAcctotET;// total ET distribution
	
	TH2F *fHistConvElectronAccEtaEET;//ET - Eta vs E   
	TH2F *fHistConvElectronAccEtaPtET;//ET - Eta vs pt  
	TH2F *fHistConvElectronAccEtaET;//  ET - Eta
	TH2F *fHistConvElectronAccEtaE;//  multiplicity - Eta vs E
	TH2F *fHistConvElectronAccEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistConvElectronAcctotET;//  total ET distribution
	
	TH2F *fHistScatElectronAccEtaEET;//ET - Eta vs E   
	TH2F *fHistScatElectronAccEtaPtET;//ET - Eta vs pt  
	TH2F *fHistScatElectronAccEtaET;//  ET - Eta
	TH2F *fHistScatElectronAccEtaE;//  multiplicity - Eta vs E
	TH2F *fHistScatElectronAccEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistScatElectronAcctotET;//  total ET distribution
	
	// *******************
	// total electron ET inside EMCal acceptance
	// *******************
	TH1F *fHistTotElectronAcctotET;//total ET distribution

	// *******************
	// gamma ET inside EMCal acceptance
	// *******************
	TH2F *fHistGammaAccEtaEET;//ET - Eta vs E   
	TH2F *fHistGammaAccEtaPtET;//ET - Eta vs pt  
	TH2F *fHistGammaAccEtaET;//  ET - Eta
	TH2F *fHistGammaAccEtaE;//  multiplicity - Eta vs E
	TH2F *fHistGammaAccEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistGammaAcctotET;//  total ET distribution
	
	TH2F *fHistAnnihGammaAccEtaEET;//ET - Eta vs E 
	TH2F *fHistAnnihGammaAccEtaPtET;//ET - Eta vs pt  
	TH2F *fHistAnnihGammaAccEtaET;//  ET - Eta
	TH2F *fHistAnnihGammaAccEtaE;//  multiplicity - Eta vs E
	TH2F *fHistAnnihGammaAccEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistAnnihGammaAcctotET;//  total ET distribution
	
	TH2F *fHistScatGammaAccEtaEET;//ET - Eta vs E  
	TH2F *fHistScatGammaAccEtaPtET;//ET - Eta vs pt  
	TH2F *fHistScatGammaAccEtaET;//  ET - Eta
	TH2F *fHistScatGammaAccEtaE;//  multiplicity - Eta vs E
	TH2F *fHistScatGammaAccEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistScatGammaAcctotET;//  total ET distribution
	
	TH2F *fHistConvGammaAccEtaEET;//ET - Eta vs E   
	TH2F *fHistConvGammaAccEtaPtET;//ET - Eta vs pt 
	TH2F *fHistConvGammaAccEtaET;//  ET - Eta
	TH2F *fHistConvGammaAccEtaE;//  multiplicity - Eta vs E
	TH2F *fHistConvGammaAccEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistConvGammaAcctotET;//  total ET distribution
	
	TH2F *fHistNonConvGammaAccEtaEET;//ET - Eta vs E   
	TH2F *fHistNonConvGammaAccEtaPtET;//ET - Eta vs pt 
	TH2F *fHistNonConvGammaAccEtaET;//  ET - Eta
	TH2F *fHistNonConvGammaAccEtaE;//  multiplicity - Eta vs E
	TH2F *fHistNonConvGammaAccEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistNonConvGammaAcctotET;//  total ET distribution
	
	// *******************
	// total gamma ET inside EMCal acceptance
	// *******************
	TH1F *fHistTotGammaAcctotET;//total ET distribution

	// *******************
	// total electromagnetic ET inside EMCal acceptance
	// *******************
	TH1F *fHistTotEMAcctotET;//total ET distribution

	// non-primary electromagnetic ET
	TH2F *fHistNPPElectronAccEtaEET;//ET - Eta vs E 
	TH2F *fHistNPPElectronAccEtaPtET;//ET - Eta vs pt 
	TH2F *fHistNPPElectronAccEtaE;// multiplicity - Eta vs E
	TH2F *fHistNPPElectronAccEtaPt;// multiplicity - Eta vs pt
	
	TH2F *fHistNPPGammaAccEtaEET;//ET - Eta vs E  
	TH2F *fHistNPPGammaAccEtaPtET;//ET - Eta vs pt 
	TH2F *fHistNPPGammaAccEtaE;// multiplicity - Eta vs E
	TH2F *fHistNPPGammaAccEtaPt;// 	multiplicity - Eta vs pt
	
	// *******************
	// electron ET reconstructed in EMCal
	// *******************
	TH2F *fHistElectronRecEtaEET;//ET - Eta vs E  
	TH2F *fHistElectronRecEtaPtET;//ET - Eta vs pt
	TH2F *fHistElectronRecEtaET;// ET - Eta
	TH2F *fHistElectronRecEtaE;// multiplicity - Eta vs E
	TH2F *fHistElectronRecEtaPt;// multiplicity - Eta vs pt
	TH1F *fHistElectronRectotET;// total ET distribution
	
	TH2F *fHistConvElectronRecEtaEET;//ET - Eta vs E   
	TH2F *fHistConvElectronRecEtaPtET;//ET - Eta vs pt 
	TH2F *fHistConvElectronRecEtaET;//  ET - Eta
	TH2F *fHistConvElectronRecEtaE;//  multiplicity - Eta vs E
	TH2F *fHistConvElectronRecEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistConvElectronRectotET;//  total ET distribution
	
	TH2F *fHistScatElectronRecEtaEET;//ET - Eta vs E   
	TH2F *fHistScatElectronRecEtaPtET;//ET - Eta vs pt  
	TH2F *fHistScatElectronRecEtaET;//  ET - Eta
	TH2F *fHistScatElectronRecEtaE;//  multiplicity - Eta vs E
	TH2F *fHistScatElectronRecEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistScatElectronRectotET;//  total ET distribution
	
	// *******************
	// total Electron ET reconstructed in EMCal
	// *******************
	TH1F *fHistTotElectronRectotET;//total ET distribution

	// *******************
	// gamma ET reconstructed in EMCal
	// *******************
	TH2F *fHistGammaRecEtaEET;//ET - Eta vs E   
	TH2F *fHistGammaRecEtaPtET;//ET - Eta vs pt  
	TH2F *fHistGammaRecEtaET;//  ET - Eta
	TH2F *fHistGammaRecEtaE;//  multiplicity - Eta vs E
	TH2F *fHistGammaRecEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistGammaRectotET;//  total ET distribution
	
	TH2F *fHistAnnihGammaRecEtaEET;//ET - Eta vs E   
	TH2F *fHistAnnihGammaRecEtaPtET;//ET - Eta vs pt  
	TH2F *fHistAnnihGammaRecEtaET;//  ET - Eta
	TH2F *fHistAnnihGammaRecEtaE;//  multiplicity - Eta vs E
	TH2F *fHistAnnihGammaRecEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistAnnihGammaRectotET;//  total ET distribution
	
	TH2F *fHistScatGammaRecEtaEET;//ET - Eta vs E   
	TH2F *fHistScatGammaRecEtaPtET;//ET - Eta vs pt  
	TH2F *fHistScatGammaRecEtaET;//  ET - Eta
	TH2F *fHistScatGammaRecEtaE;//  multiplicity - Eta vs E
	TH2F *fHistScatGammaRecEtaPt;//  multiplicity - Eta vs pt
	TH1F *fHistScatGammaRectotET;//  total ET distribution

	// *******************
	// total gamma ET reconstructed in EMCal
	// *******************
	TH1F *fHistTotGammaRectotET;//total ET distribution

	// *******************
	// total EM ET reconstructed in EMCal
	// *******************
	TH1F *fHistTotEMRectotET;//total ET distribution

	// non-primary electromagnetic ET
	TH2F *fHistNPPElectronRecEtaEET;//ET - Eta vs E  
	TH2F *fHistNPPElectronRecEtaPtET;//ET - Eta vs pt 
	TH2F *fHistNPPElectronRecEtaET;// ET - Eta
	TH2F *fHistNPPElectronRecEtaE;// multiplicity - Eta vs E
	TH2F *fHistNPPElectronRecEtaPt;// multiplicity - Eta vs pt
	TH1F *fHistNPPElectronRectotET;// total ET distribution
	
	TH2F *fHistNPPGammaRecEtaEET;//ET - Eta vs E  
	TH2F *fHistNPPGammaRecEtaPtET;//ET - Eta vs pt 
	TH2F *fHistNPPGammaRecEtaET;// ET - Eta
	TH2F *fHistNPPGammaRecEtaE;// multiplicity - Eta vs E
	TH2F *fHistNPPGammaRecEtaPt;// multiplicity - Eta vs pt
	TH1F *fHistNPPGammaRectotET;// total ET distribution
	
	TH1F *fHistTotNPPEMRectotET;//total ET distribution

	TH2F *fHistNPPPi0GammaRecEtaEET;//ET - Eta vs E  
	TH2F *fHistNPPPi0GammaRecEtaPtET;//ET - Eta vs pt 
	TH2F *fHistNPPPi0GammaRecEtaET;// ET - Eta
	TH2F *fHistNPPPi0GammaRecEtaE;// multiplicity - Eta vs E
	TH2F *fHistNPPPi0GammaRecEtaPt;// multiplicity - Eta vs pt
	TH1F *fHistNPPPi0GammaRectotET;// total ET distribution
	
	// *******************
	// muon ET (+ and -)
	// *******************
	TH2F *fHistMuonEtaEET;//ET - Eta vs E  
	TH2F *fHistMuonAccEtaEET;//ET - Eta vs E 
	TH2F *fHistMuonRecEtaEET;//ET - Eta vs E 
	TH2F *fHistMuonMatchEtaEET;//ET - Eta vs E 

	TH2F *fHistMuonEtaPtET;// ET - Eta vs pt
	TH2F *fHistMuonAccEtaPtET;// ET - Eta vs pt
	TH2F *fHistMuonRecEtaPtET;// ET - Eta vs pt
	TH2F *fHistMuonMatchEtaPtET;// ET - Eta vs pt

	TH2F *fHistMuonEtaET;// ET - Eta
	TH2F *fHistMuonAccEtaET;// ET - Eta
	TH2F *fHistMuonRecEtaET;// ET - Eta
	TH2F *fHistMuonMatchEtaET;// ET - Eta
	
	TH2F *fHistMuonEtaE;// multiplicity - Eta vs E
	TH2F *fHistMuonAccEtaE;// multiplicity - Eta vs E
	TH2F *fHistMuonRecEtaE;// multiplicity - Eta vs E
	TH2F *fHistMuonMatchEtaE;// multiplicity - Eta vs E
	
	TH2F *fHistMuonEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistMuonAccEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistMuonRecEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistMuonMatchEtaPt;// multiplicity - Eta vs pt
	
	TH1F *fHistMuontotET;// total ET distribution
	TH1F *fHistMuonAcctotET;// total ET distribution
	TH1F *fHistMuonRectotET;// total ET distribution
	TH1F *fHistMuonMatchtotET;// total ET distribution
	
	TH1F *fHistMuonRectotETDep;//total deposited ET distribution
	TH1F *fHistMuonMatchtotETDep;// total deposited ET distribution
	
	TH2F *fHistMuonRecEtaEDepETDep;// ET deposited - Eta vs E deposited
	TH2F *fHistMuonMatchEtaEDepETDep;// ET deposited - Eta vs E deposited

	TH2F *fHistMuonRecEtaPtETDep;// ET deposited - Eta vs pt
	TH2F *fHistMuonMatchEtaPtETDep;// ET deposited - Eta vs pt
	
	TH2F *fHistMuonRecEtaETDep;// ET deposited - Eta
	TH2F *fHistMuonMatchEtaETDep;// ET deposited - Eta

	TH2F *fHistMuonRecResEET;// ET - track matching residual vs E
	TH2F *fHistMuonRecResPtET;// ET - track matching residual vs pt
	TH2F *fHistMuonRecResE;// multiplicity - track matching residual vs E
	TH2F *fHistMuonRecResPt;// multiplicity - track matching residual vs pt 
	TH2F *fHistMuonRecResEDepETDep;// ET deposited - track matching residual vs E deposited
	TH2F *fHistMuonRecResPtETDep;// ET deposited - track matching residual vs pt
	
	// *******************
	// pion ET (+ and -)
	// *******************
	TH2F *fHistPionEtaEET;//ET - Eta vs E  
	TH2F *fHistPionAccEtaEET;//ET - Eta vs E 
	TH2F *fHistPionRecEtaEET;//ET - Eta vs E 
	TH2F *fHistPionMatchEtaEET;//ET - Eta vs E 
	
	TH2F *fHistPionEtaPtET;// ET - Eta vs pt
	TH2F *fHistPionAccEtaPtET;// ET - Eta vs pt
	TH2F *fHistPionRecEtaPtET;// ET - Eta vs pt
	TH2F *fHistPionMatchEtaPtET;// ET - Eta vs pt
	
	TH2F *fHistPionEtaET;// ET - Eta
	TH2F *fHistPionAccEtaET;// ET - Eta
	TH2F *fHistPionRecEtaET;// ET - Eta
	TH2F *fHistPionMatchEtaET;// ET - Eta
	
	TH2F *fHistPionEtaE;// multiplicity - Eta vs E
	TH2F *fHistPionAccEtaE;// multiplicity - Eta vs E
	TH2F *fHistPionRecEtaE;// multiplicity - Eta vs E
	TH2F *fHistPionMatchEtaE;// multiplicity - Eta vs E
	
	TH2F *fHistPionEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistPionAccEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistPionRecEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistPionMatchEtaPt;// multiplicity - Eta vs pt
	
	TH1F *fHistPiontotET;// total ET distribution
	TH1F *fHistPionAcctotET;// total ET distribution
	TH1F *fHistPionRectotET;// total ET distribution
	TH1F *fHistPionMatchtotET;// total ET distribution
	
	TH1F *fHistPionRectotETDep;// total deposited ET distribution
	TH1F *fHistPionMatchtotETDep;// total deposited ET distribution
	
	TH2F *fHistPionRecEtaEDepETDep;// ET deposited - Eta vs E deposited
	TH2F *fHistPionMatchEtaEDepETDep;// ET deposited - Eta vs E deposited

	TH2F *fHistPionRecEtaPtETDep;// ET deposited - Eta vs pt
	TH2F *fHistPionMatchEtaPtETDep;// ET deposited - Eta vs pt
	
	TH2F *fHistPionRecEtaETDep;// ET deposited - Eta
	TH2F *fHistPionMatchEtaETDep;// ET deposited - Eta
	
	TH2F *fHistPionRecResEET;// ET - track matching residual vs E
	TH2F *fHistPionRecResPtET;// ET - track matching residual vs pt
	TH2F *fHistPionRecResE;// multiplicity - track matching residual vs E
	TH2F *fHistPionRecResPt;// multiplicity - track matching residual vs pt
	TH2F *fHistPionRecResEDepETDep;// ET deposited - track matching residual vs E deposited
	TH2F *fHistPionRecResPtETDep;// ET deposited - track matching residual vs pt
	
	// *******************
	// charged kaon (+ and -) ET
	// *******************
	TH2F *fHistKaonEtaEET;//ET - Eta vs E 
	TH2F *fHistKaonAccEtaEET;//ET - Eta vs E 
	TH2F *fHistKaonRecEtaEET;//ET - Eta vs E 
	TH2F *fHistKaonMatchEtaEET;//ET - Eta vs E 
	
	TH2F *fHistKaonEtaPtET;// ET - Eta vs pt
	TH2F *fHistKaonAccEtaPtET;// ET - Eta vs pt
	TH2F *fHistKaonRecEtaPtET;// ET - Eta vs pt
	TH2F *fHistKaonMatchEtaPtET;// ET - Eta vs pt
	
	TH2F *fHistKaonEtaET;// ET - Eta
	TH2F *fHistKaonAccEtaET;// ET - Eta
	TH2F *fHistKaonRecEtaET;// ET - Eta
	TH2F *fHistKaonMatchEtaET;// ET - Eta
	
	TH2F *fHistKaonEtaE;// multiplicity - Eta vs E
	TH2F *fHistKaonAccEtaE;// multiplicity - Eta vs E
	TH2F *fHistKaonRecEtaE;// multiplicity - Eta vs E
	TH2F *fHistKaonMatchEtaE;// multiplicity - Eta vs E
	
	TH2F *fHistKaonEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistKaonAccEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistKaonRecEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistKaonMatchEtaPt;// multiplicity - Eta vs pt

	TH1F *fHistKaontotET;// total ET distribution
	TH1F *fHistKaonAcctotET;// total ET distribution
	TH1F *fHistKaonRectotET;// total ET distribution
	TH1F *fHistKaonMatchtotET;// total ET distribution
	
	TH1F *fHistKaonRectotETDep;// total deposited ET distribution
	TH1F *fHistKaonMatchtotETDep;// total deposited ET distribution
	
	TH2F *fHistKaonRecEtaEDepETDep;// ET deposited - Eta vs E deposited
	TH2F *fHistKaonMatchEtaEDepETDep;// ET deposited - Eta vs E deposited

	TH2F *fHistKaonRecEtaPtETDep;// ET deposited - Eta vs pt
	TH2F *fHistKaonMatchEtaPtETDep;// ET deposited - Eta vs pt
	
	TH2F *fHistKaonRecEtaETDep;// ET deposited - Eta
	TH2F *fHistKaonMatchEtaETDep;// ET deposited - Eta
	
	TH2F *fHistKaonRecResEET;// ET - track matching residual vs E
	TH2F *fHistKaonRecResPtET;// ET - track matching residual vs pt
	TH2F *fHistKaonRecResE;// multiplicity - track matching residual vs E
	TH2F *fHistKaonRecResPt;// multiplicity - track matching residual vs pt
	TH2F *fHistKaonRecResEDepETDep;// ET deposited - track matching residual vs E deposited
	TH2F *fHistKaonRecResPtETDep;// ET deposited - track matching residual vs pt
	
	// *******************
	// proton (anti) ET
	// *******************
	TH2F *fHistProtonEtaEET;//ET - Eta vs E 
	TH2F *fHistProtonAccEtaEET;//ET - Eta vs E 
	TH2F *fHistProtonRecEtaEET;//ET - Eta vs E 
	TH2F *fHistProtonMatchEtaEET;//ET - Eta vs E 
	
	TH2F *fHistProtonEtaPtET;// ET - Eta vs pt
	TH2F *fHistProtonAccEtaPtET;// ET - Eta vs pt
	TH2F *fHistProtonRecEtaPtET;// ET - Eta vs pt
	TH2F *fHistProtonMatchEtaPtET;// ET - Eta vs pt
	
	TH2F *fHistProtonEtaET;// ET - Eta
	TH2F *fHistProtonAccEtaET;// ET - Eta
	TH2F *fHistProtonRecEtaET;// ET - Eta
	TH2F *fHistProtonMatchEtaET;// ET - Eta
	
	TH2F *fHistProtonEtaE;// multiplicity - Eta vs E
	TH2F *fHistProtonAccEtaE;// multiplicity - Eta vs E
	TH2F *fHistProtonRecEtaE;// multiplicity - Eta vs E
	TH2F *fHistProtonMatchEtaE;// multiplicity - Eta vs E
	
	TH2F *fHistProtonEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistProtonAccEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistProtonRecEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistProtonMatchEtaPt;// multiplicity - Eta vs pt

	TH1F *fHistProtontotET;// total ET distribution
	TH1F *fHistProtonAcctotET;// total ET distribution
	TH1F *fHistProtonRectotET;// total ET distribution
	TH1F *fHistProtonMatchtotET;// total ET distribution
	
	TH1F *fHistProtonRectotETDep;// total deposited ET distribution
	TH1F *fHistProtonMatchtotETDep;// total deposited ET distribution
	
	TH2F *fHistProtonRecEtaEDepETDep;// ET deposited - Eta vs E deposited
	TH2F *fHistProtonMatchEtaEDepETDep;// ET deposited - Eta vs E deposited
	
	TH2F *fHistProtonRecEtaPtETDep;// ET deposited - Eta vs pt
	TH2F *fHistProtonMatchEtaPtETDep;// ET deposited - Eta vs pt
	
	TH2F *fHistProtonRecEtaETDep;// ET deposited - Eta
	TH2F *fHistProtonMatchEtaETDep;// ET deposited - Eta

	TH2F *fHistProtonRecResEET;// ET - track matching residual vs E
	TH2F *fHistProtonRecResPtET;// ET - track matching residual vs pt
	TH2F *fHistProtonRecResE;// multiplicity - track matching residual vs E
	TH2F *fHistProtonRecResPt;// multiplicity - track matching residual vs pt
	TH2F *fHistProtonRecResEDepETDep;// ET deposited - track matching residual vs E deposited
	TH2F *fHistProtonRecResPtETDep;// ET deposited - track matching residual vs pt
	
	// *******************
	// total charged ET
	// *******************
	TH1F *fHistTotChargedtotET;//total ET distribution
	TH1F *fHistTotChargedAcctotET;//total ET distribution
	TH1F *fHistTotChargedRectotET;//total ET distribution
	TH1F *fHistTotChargedRectotETDep;//total deposited ET distribution
	TH1F *fHistTotChargedMatchtotET;//total ET distribution
	TH1F *fHistTotChargedMatchtotETDep;//total deposited ET distribution
	
	// *******************
	// neutron (anti) ET
	// *******************
	TH2F *fHistNeutronEtaEET;//ET - Eta vs E 
	TH2F *fHistNeutronAccEtaEET;//ET - Eta vs E 
	TH2F *fHistNeutronRecEtaEET;//ET - Eta vs E 
	
	TH2F *fHistNeutronEtaPtET;// ET - Eta vs pt
	TH2F *fHistNeutronAccEtaPtET;// ET - Eta vs pt
	TH2F *fHistNeutronRecEtaPtET;// ET - Eta vs pt
	
	TH2F *fHistNeutronEtaET;// ET - Eta
	TH2F *fHistNeutronAccEtaET;// ET - Eta
	TH2F *fHistNeutronRecEtaET;// ET - Eta
	
	TH2F *fHistNeutronEtaE;// multiplicity - Eta vs E
	TH2F *fHistNeutronAccEtaE;// multiplicity - Eta vs E
	TH2F *fHistNeutronRecEtaE;// multiplicity - Eta vs E
	
	TH2F *fHistNeutronEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistNeutronAccEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistNeutronRecEtaPt;// multiplicity - Eta vs pt
	
	TH1F *fHistNeutrontotET;// total ET distribution
	TH1F *fHistNeutronAcctotET;// total ET distribution
	TH1F *fHistNeutronRectotET;// total ET distribution
	TH1F *fHistNeutronRectotETDep;// total deposited ET distribution
	
	TH2F *fHistNeutronRecEtaEDepETDep;// ET deposited - Eta vs E deposited
	TH2F *fHistNeutronRecEtaETDep;// ET deposited - Eta
	
	TH2F *fHistNeutronRecEtaPtETDep;// ET deposited - Eta vs pt
		
	// *******************
	// neutral kaon ET
	// *******************
	TH2F *fHistK0EtaEET;//ET - Eta vs E 
	TH2F *fHistK0RecEtaEET;//ET - Eta vs E 
	
	TH2F *fHistK0EtaPtET;// ET - Eta vs pt
	TH2F *fHistK0RecEtaPtET;// ET - Eta vs pt
	
	TH2F *fHistK0EtaET;// ET - Eta
	TH2F *fHistK0RecEtaET;// ET - Eta
	
	TH2F *fHistK0EtaE;// multiplicity - Eta vs E
	TH2F *fHistK0RecEtaE;// multiplicity - Eta vs E
	
	TH2F *fHistK0EtaPt;// multiplicity - Eta vs pt
	TH2F *fHistK0RecEtaPt;// multiplicity - Eta vs pt

	TH1F *fHistK0totET;// total ET distribution
	TH1F *fHistK0RectotET;// total ET distribution
	
	TH1F *fHistK0RectotETDep;// total deposited ET distribution
	
	TH2F *fHistK0RecEtaEDepETDep;// ET deposited - Eta vs E deposited
	TH2F *fHistK0RecEtaETDep;// ET deposited - Eta
	
	TH2F *fHistK0RecEtaPtETDep;// ET deposited - Eta vs pt
		
	// *******************
	// Lambda(anti) ET
	// *******************
	TH2F *fHistLambdaEtaEET;//ET - Eta vs E 
	TH2F *fHistLambdaRecEtaEET;//ET - Eta vs E 
	
	TH2F *fHistLambdaEtaPtET;// ET - Eta vs pt
	TH2F *fHistLambdaRecEtaPtET;// ET - Eta vs pt
	
	TH2F *fHistLambdaEtaET;// ET - Eta
	TH2F *fHistLambdaRecEtaET;// ET - Eta
	
	TH2F *fHistLambdaEtaE;// multiplicity - Eta vs E
	TH2F *fHistLambdaRecEtaE;// multiplicity - Eta vs E
	
	TH2F *fHistLambdaEtaPt;// multiplicity - Eta vs pt
	TH2F *fHistLambdaRecEtaPt;// multiplicity - Eta vs pt
	
	TH1F *fHistLambdatotET;// total ET distribution
	TH1F *fHistLambdaRectotET;// total ET distribution
	
	TH1F *fHistLambdaRectotETDep;// total deposited ET distribution
	
	TH2F *fHistLambdaRecEtaEDepETDep;// ET deposited - Eta vs E deposited
	TH2F *fHistLambdaRecEtaETDep;// ET deposited - Eta
	
	TH2F *fHistLambdaRecEtaPtETDep;// ET deposited - Eta vs pt

	// *******************
	// total neutral ET
	// *******************
	TH1F *fHistTotNeutraltotET;//total ET distribution
	TH1F *fHistTotNeutralRectotET;//total ET distribution
	TH1F *fHistTotNeutralRectotETDep;//total deposited ET distribution
	
	// *******************
	// total ET
	// *******************
	TH1F *fHistTotaltotET;//total ET distribution
	TH1F *fHistTotalAcctotET;//total ET distribution
	TH1F *fHistTotalRectotET;//total ET distribution
	TH1F *fHistTotalRectotETDep;//total deposited ET distribution
	
	// *******************
	// some checks
	// *******************

	// check produced electrons
	TH1F *fHistElectronFirstMother;// first mother ID
	TH2F *fHistElectronFirstMotherXY;// first mother XY position 
	TH1F *fHistElectronNDaughters;// number of daughters 
	TH1F *fHistElectronDaughters;// daughters ID
	TH2F *fHistElectronDaughtersXY;// daughters XY position

	TH1F *fHistElectronFirstMotherAcc;// first mother ID
	TH2F *fHistElectronFirstMotherXYAcc;// first mother XY position 
	TH1F *fHistElectronNDaughtersAcc;// number of daughters
	TH1F *fHistElectronDaughtersAcc;// daughters ID
	TH2F *fHistElectronDaughtersXYAcc;// daughters XY position

	TH1F *fHistElectronFirstMotherRec;// first mother ID 
	TH2F *fHistElectronFirstMotherXYRec;//  first mother XY position
	TH1F *fHistElectronNDaughtersRec;// number of daughters
	TH1F *fHistElectronDaughtersRec;// daughters ID
	TH2F *fHistElectronDaughtersXYRec;// daughters XY position

	TH1F *fHistNPPElectronFirstMother;//  first mother ID 
	TH2F *fHistNPPElectronFirstMotherXY;// first mother XY position
	TH1F *fHistNPPElectronNDaughters;// number of daughters
	TH1F *fHistNPPElectronDaughters;// daughters ID
	TH2F *fHistNPPElectronDaughtersXY;// daughters XY position
	
	TH1F *fHistNPPElectronFirstMotherAcc;//  first mother ID 
	TH2F *fHistNPPElectronFirstMotherXYAcc;//  first mother XY position
	TH1F *fHistNPPElectronNDaughtersAcc;// number of daughters
	TH1F *fHistNPPElectronDaughtersAcc;// daughters ID
	TH2F *fHistNPPElectronDaughtersXYAcc;// daughters XY position
	
	TH1F *fHistNPPElectronFirstMotherRec;// first mother ID  
	TH2F *fHistNPPElectronFirstMotherXYRec;//  first mother XY position
	TH1F *fHistNPPElectronNDaughtersRec;// number of daughters
	TH1F *fHistNPPElectronDaughtersRec;// daughters ID
	TH2F *fHistNPPElectronDaughtersXYRec;// daughters XY position
	
	// check produced gammas
	TH1F *fHistGammaFirstMother;// first mother ID 
	TH2F *fHistGammaFirstMotherXY;// first mother XY position
	TH1F *fHistGammaNDaughters;// number of daughters
	TH1F *fHistGammaDaughters;// daughters ID
	TH2F *fHistGammaDaughtersXY;// daughters XY position
	TH2F *fHistConvGammaDaughtersXY;// daughters XY position
	TH2F *fHistNonConvGammaDaughtersXY;// daughters XY position
	
	TH1F *fHistGammaFirstMotherAcc;// first mother ID  
	TH2F *fHistGammaFirstMotherXYAcc;//  first mother XY position
	TH1F *fHistGammaNDaughtersAcc;// number of daughters
	TH1F *fHistGammaDaughtersAcc;// daughters ID
	TH2F *fHistGammaDaughtersXYAcc;// daughters XY position
	TH2F *fHistConvGammaDaughtersXYAcc;// daughters XY position
	TH2F *fHistNonConvGammaDaughtersXYAcc;// daughters XY position
	
	TH1F *fHistGammaFirstMotherRec;// first mother ID  
	TH2F *fHistGammaFirstMotherXYRec;//  first mother XY position
	TH1F *fHistGammaNDaughtersRec;// number of daughters
	TH1F *fHistGammaDaughtersRec;// daughters ID
	TH2F *fHistGammaDaughtersXYRec;// daughters XY position
	TH2F *fHistConvGammaDaughtersXYRec;// daughters XY position
	TH2F *fHistNonConvGammaDaughtersXYRec;// daughters XY position
	
	TH1F *fHistNPPGammaFirstMother;// first mother ID 
	TH2F *fHistNPPGammaFirstMotherXY;// first mother XY position
	TH1F *fHistNPPGammaNDaughters;// number of daughters
	TH1F *fHistNPPGammaDaughters;// daughters ID
	TH2F *fHistNPPGammaDaughtersXY;// daughters XY position
	
	TH1F *fHistNPPGammaFirstMotherAcc;// first mother ID  
	TH2F *fHistNPPGammaFirstMotherXYAcc;//  first mother XY position
	TH1F *fHistNPPGammaNDaughtersAcc;// number of daughters
	TH1F *fHistNPPGammaDaughtersAcc;// daughters ID
	TH2F *fHistNPPGammaDaughtersXYAcc;// daughters XY position
	
	TH1F *fHistNPPGammaFirstMotherRec;// first mother ID  
	TH2F *fHistNPPGammaFirstMotherXYRec;//  first mother XY position
	TH1F *fHistNPPGammaNDaughtersRec;// number of daughters
	TH1F *fHistNPPGammaDaughtersRec;// daughters ID
	TH2F *fHistNPPGammaDaughtersXYRec;// daughters XY position

	//check projections
	TH2F *fHistAllERecEMC;// E reconstructed vs E MC	
	TH2F *fHistAllPtRecPtMC;// pt reconstructed vs pt MC
	TH2F *fHistElectronERecEMC;// E reconstructed vs E MC	
	TH2F *fHistGammaERecEMC;// E reconstructed vs E MC
	
	TH2F *fHistChargedRes;// charged particle track matching residual
	TH2F *fHistChargedRes2;// charged particle track matching residual
	TH2F *fHistChargedRes3;// charged particle track matching residual
	TH2F *fHistNeutralRes;// neutral particle track matching residual
	TH2F *fHistElectronRes;// electron track matching residual
	TH2F *fHistGammaRes;// gamma track matching residual
	
	TH2F *fHistIsInAcc;// EMCal acceptance check
	
 private:

  //Declare it private to avoid compilation warning
    AliAnalysisEmEtMonteCarlo & operator = (const AliAnalysisEmEtMonteCarlo & g) ;//cpy assignment
    AliAnalysisEmEtMonteCarlo(const AliAnalysisEmEtMonteCarlo & g) ; // cpy ctor
    ClassDef(AliAnalysisEmEtMonteCarlo, 1);
};

#endif //ALIANALYSISEMETMONTECARLO_H
