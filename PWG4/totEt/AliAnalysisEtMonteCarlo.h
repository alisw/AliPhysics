#ifndef ALIANALYSISETMONTECARLO_H
#define ALIANALYSISETMONTECARLO_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis
//  - MC output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEt.h"
class TParticle;
class TH3F;
//class AliMCEvent;
//class AliESDEvent;

class AliAnalysisEtMonteCarlo : public AliAnalysisEt
{

public:
   
  AliAnalysisEtMonteCarlo();
  virtual ~AliAnalysisEtMonteCarlo();

    virtual Int_t AnalyseEvent(AliVEvent* event);
	virtual Int_t AnalyseEvent(AliVEvent* event, AliVEvent* event2);
	//virtual Int_t AnalyseEvent(AliMCEvent* event, AliESDEvent* event2);

    virtual void Init();
    virtual void ResetEventValues();
    virtual void CreateHistograms();
    virtual void FillOutputList(TList* list);

    virtual void FillHistograms();

protected:

    virtual bool TrackHitsCalorimeter(TParticle *part, Double_t magField=0.5);

protected:

    Double_t fImpactParameter; // b(fm), for Hijing; 0 otherwise
    Int_t fNcoll; // Ncoll, for Hijing; 1 otherwise
    Int_t fNpart; // Ncoll, for Hijing; 2 otherwise

	TH2F *fHistPrimElectronEtaEET; // electrons that are physical primaries 
	TH2F *fHistSecElectronEtaEET; // electrons that are NOT physical primaries and DO NOT come from gammas (that are physical primaries) conversion 
	TH2F *fHistConvElectronEtaEET; // electrons that are NOT physical primaries and come from gammas (that are physical primaries) conversion 
	TH2F *fHistPrimGammaEtaEET; // gammas (physical primaries) that DO NOT come from pi0, eta or omega0 decays 
	TH2F *fHistPion0GammaEtaEET; // gammas (physical primaries) that come from pi0 decays 
	TH2F *fHistEtaGammaEtaEET; // gammas (physical primaries) that come from eta decays 
	TH2F *fHistOmega0GammaEtaEET; // gammas (physical primaries) that come from omega0 decays 
	TH2F *fHistSecGammaEtaEET; // gammas that are not physical primaries 
	
	TH2F *fHistPrimElectronEtaE; // simple couting of MC input
	TH2F *fHistSecElectronEtaE;//Histogram of secondary elections
	TH2F *fHistConvElectronEtaE;//conversion electrons
	TH2F *fHistPrimGammaEtaE; //primary gammas
	TH2F *fHistPion0GammaEtaE;//gammas from pions
	TH2F *fHistEtaGammaEtaE;//gammas from etas
	TH2F *fHistOmega0GammaEtaE;//gammas from omegas
	TH2F *fHistSecGammaEtaE; //secondary gammas
	
	TH2F *fHistPrimElectronEtaERec; // simple couting of recosntructed MC input
	TH2F *fHistSecElectronEtaERec;//secondary electrons
	TH2F *fHistConvElectronEtaERec;//conversion electrons
	TH2F *fHistPrimGammaEtaERec; //primary gammas
	TH2F *fHistSecGammaEtaERec; //secondary gammas
	TH2F *fHistPion0GammaEtaERec;//pion gammas
	TH2F *fHistEtaGammaEtaERec;//eta gammas
	TH2F *fHistOmega0GammaEtaERec;//omega gammas
	
	TH2F *fHistAllERecEMC;	// compare energy directly from MC and the recosntructed one 
	TH2F *fHistGammaERecEMC;//gammas 	
	TH2F *fHistElectronERecEMC;	//gammas
	
	TH1F *fHistElectronFirstMother; // distribution of PDG code of electron's (physical primary) first mothers 
	TH1F *fHistElectronLastMother; // distribution of PDG code of electron's (physical primary) last (or second) mothers 
	TH1F *fHistElectronFirstMotherEtaAcc; // distribution of PDG code of electron's (physical primary and inside acceptance) first mothers 
	TH1F *fHistElectronFirstMotherNPP; // distribution of PDG code of electron's (NON physical primary) first mothers 
	TH1F *fHistElectronFirstMotherNPPAcc; // distribution of PDG code of electron's (NON physical primary and inside acceptance) first mothers 

	TH1F *fHistGammaFirstMother; // same as above but for gammas 
	TH1F *fHistGammaLastMother;//enter comment here
	TH1F *fHistGammaFirstMotherEtaAcc;//enter comment here
	TH1F *fHistGammaFirstMotherNPP;//enter comment here
	TH1F *fHistGammaFirstMotherNPPAcc;//enter comment here
	
	TH3F *fHistDecayVertexNonRemovedCharged; // Decay vertex for non-removed charged particles
	TH3F *fHistDecayVertexRemovedCharged; // Decay vertex for non-removed charged particles
	TH3F *fHistDecayVertexNonRemovedNeutral; // Decay vertex for non-removed charged particles
	TH3F *fHistDecayVertexRemovedNeutral; // Decay vertex for non-removed charged particles
	
	TH2F *fHistRemovedOrNot; // If charged/neutral particles were removed or not
	
	TH2F *fHistEtNonRemovedProtons; 
	TH2F *fHistEtNonRemovedAntiProtons; 
	TH2F *fHistEtNonRemovedPiPlus; 
	TH2F *fHistEtNonRemovedPiMinus; 
	TH2F *fHistEtNonRemovedKaonPlus; 
	TH2F *fHistEtNonRemovedKaonMinus; 
	TH2F *fHistEtNonRemovedK0s; 
	TH2F *fHistEtNonRemovedLambdas; 
	TH2F *fHistEtNonRemovedElectrons; 
	TH2F *fHistEtNonRemovedPositrons; 
	TH2F *fHistEtNonRemovedMuPlus;
	TH2F *fHistEtNonRemovedMuMinus;
	TH2F *fHistEtNonRemovedNeutrons;
	TH2F *fHistEtNonRemovedAntiNeutrons;
	TH2F *fHistEtNonRemovedGammas;
	TH2F *fHistEtNonRemovedGammasFromPi0;
	
	TH2F *fHistEtRemovedGammas;
	TH2F *fHistEtRemovedNeutrons;
	TH2F *fHistEtRemovedAntiNeutrons;
	
	
	TH2F *fHistMultNonRemovedProtons; 
	TH2F *fHistMultNonRemovedAntiProtons; 
	TH2F *fHistMultNonRemovedPiPlus; 
	TH2F *fHistMultNonRemovedPiMinus; 
	TH2F *fHistMultNonRemovedKaonPlus; 
	TH2F *fHistMultNonRemovedKaonMinus; 
	TH2F *fHistMultNonRemovedK0s; 
	TH2F *fHistMultNonRemovedLambdas; 
	TH2F *fHistMultNonRemovedElectrons; 
	TH2F *fHistMultNonRemovedPositrons; 
	TH2F *fHistMultNonRemovedMuPlus;
	TH2F *fHistMultNonRemovedMuMinus;
	TH2F *fHistMultNonRemovedNeutrons;
	TH2F *fHistMultNonRemovedAntiNeutrons;
	TH2F *fHistMultNonRemovedGammas;
	
	TH2F *fHistMultRemovedGammas;
	TH2F *fHistMultRemovedNeutrons;
	TH2F *fHistMultRemovedAntiNeutrons;
	
	TH2F *fHistTrackMultvsNonRemovedCharged;
	TH2F *fHistTrackMultvsNonRemovedNeutral;
	TH2F *fHistTrackMultvsRemovedGamma;
	
	TH2F *fHistClusterMultvsNonRemovedCharged;
	TH2F *fHistClusterMultvsNonRemovedNeutral;
	TH2F *fHistClusterMultvsRemovedGamma;
	
	TH2F *fHistMultvsNonRemovedChargedE;
	TH2F *fHistMultvsNonRemovedNeutralE;
	TH2F *fHistMultvsRemovedGammaE;
	
	Float_t fEtNonRemovedProtons; 
	Float_t fEtNonRemovedAntiProtons; 
	Float_t fEtNonRemovedPiPlus; 
	Float_t fEtNonRemovedPiMinus; 
	Float_t fEtNonRemovedKaonPlus; 
	Float_t fEtNonRemovedKaonMinus; 
	Float_t fEtNonRemovedK0s; 
	Float_t fEtNonRemovedLambdas; 
	Float_t fEtNonRemovedElectrons; 
	Float_t fEtNonRemovedPositrons; 
	Float_t fEtNonRemovedMuMinus;
	Float_t fEtNonRemovedMuPlus;
	Float_t fEtNonRemovedGammas;
	Float_t fEtNonRemovedGammasFromPi0;
	Float_t fEtNonRemovedNeutrons;
	Float_t fEtNonRemovedAntiNeutrons;
	
	Float_t fEtRemovedGammas;
	Float_t fEtRemovedNeutrons;
	Float_t fEtRemovedAntiNeutrons;
		
	Int_t fMultNonRemovedProtons; 
	Int_t fMultNonRemovedAntiProtons; 
	Int_t fMultNonRemovedPiPlus; 
	Int_t fMultNonRemovedPiMinus; 
	Int_t fMultNonRemovedKaonPlus; 
	Int_t fMultNonRemovedKaonMinus; 
	Int_t fMultNonRemovedK0s; 
	Int_t fMultNonRemovedLambdas; 
	Int_t fMultNonRemovedElectrons; 
	Int_t fMultNonRemovedPositrons; 
	Int_t fMultNonRemovedMuMinus;
	Int_t fMultNonRemovedMuPlus;
	Int_t fMultNonRemovedGammas;
	Int_t fMultNonRemovedNeutrons;
	Int_t fMultNonRemovedAntiNeutrons;
	
	Int_t fMultRemovedGammas;
	Int_t fMultRemovedNeutrons;
	Int_t fMultRemovedAntiNeutrons;
	
	Int_t fTrackMultInAcc;
	
	
	TH2F *fHistDxDzNonRemovedCharged;
	TH2F *fHistDxDzRemovedCharged;
	TH2F *fHistDxDzNonRemovedNeutral;
	TH2F *fHistDxDzRemovedNeutral;
	
	TH1F *fHistPiPlusMult;
	TH1F *fHistPiMinusMult;
	TH1F *fHistPiZeroMult;

	TH1F *fHistPiPlusMultAcc;
	TH1F *fHistPiMinusMultAcc;
	TH1F *fHistPiZeroMultAcc;
	
	Int_t fPiPlusMult;
	Int_t fPiMinusMult;
	Int_t fPiZeroMult;

	Int_t fPiPlusMultAcc;
	Int_t fPiMinusMultAcc;
	Int_t fPiZeroMultAcc;
	
	
	Int_t fNeutralRemoved; // number of neutral particles that where removed by track matching
	Int_t fChargedRemoved; // number of charged particles that where removed by track matching
	Int_t fChargedNotRemoved; // number of charged particles that were not removed
	Int_t fNeutralNotRemoved; // number of neutral particles that were not removed
	
	Double_t fEnergyNeutralRemoved; // energy of neutral particles that where removed by track matching
	Double_t fEnergyChargedRemoved; // energy of charged particles that where removed by track matching
	Double_t fEnergyChargedNotRemoved; // energy of charged particles that were not removed
	Double_t fEnergyNeutralNotRemoved; // energy of neutral particles that were not removed
	
	
 private:

    //Declare it private to avoid compilation warning
    AliAnalysisEtMonteCarlo & operator = (const AliAnalysisEtMonteCarlo & g) ;//cpy assignment
    AliAnalysisEtMonteCarlo(const AliAnalysisEtMonteCarlo & g) ; // cpy ctor
    ClassDef(AliAnalysisEtMonteCarlo, 2);
};

#endif // ALIANALYSISETMONTECARLO_H
