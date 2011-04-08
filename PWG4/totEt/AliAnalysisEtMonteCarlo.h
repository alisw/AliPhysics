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
	TH2F *fHistElectronERecEMC;	//electrons
	
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
	
	TH2F *fHistEtNonRemovedProtons; // enter comment here
	TH2F *fHistEtNonRemovedAntiProtons; // enter comment here
	TH2F *fHistEtNonRemovedPiPlus; // enter comment here
	TH2F *fHistEtNonRemovedPiMinus; // enter comment here
	TH2F *fHistEtNonRemovedKaonPlus; // enter comment here
	TH2F *fHistEtNonRemovedKaonMinus; // enter comment here
	TH2F *fHistEtNonRemovedK0s; // enter comment here
	TH2F *fHistEtNonRemovedLambdas; // enter comment here
	TH2F *fHistEtNonRemovedElectrons; // enter comment here
	TH2F *fHistEtNonRemovedPositrons; // enter comment here
	TH2F *fHistEtNonRemovedMuPlus; // enter comment here
	TH2F *fHistEtNonRemovedMuMinus; // enter comment here
	TH2F *fHistEtNonRemovedNeutrons; // enter comment here
	TH2F *fHistEtNonRemovedAntiNeutrons; // enter comment here
	TH2F *fHistEtNonRemovedGammas; // enter comment here
	TH2F *fHistEtNonRemovedGammasFromPi0; // enter comment here
	
	TH2F *fHistEtRemovedGammas; // enter comment here
	TH2F *fHistEtRemovedNeutrons; // enter comment here
	TH2F *fHistEtRemovedAntiNeutrons; // enter comment here
	
	
	TH2F *fHistMultNonRemovedProtons; // enter comment here 
	TH2F *fHistMultNonRemovedAntiProtons; // enter comment here 
	TH2F *fHistMultNonRemovedPiPlus; // enter comment here 
	TH2F *fHistMultNonRemovedPiMinus; // enter comment here 
	TH2F *fHistMultNonRemovedKaonPlus; // enter comment here 
	TH2F *fHistMultNonRemovedKaonMinus; // enter comment here 
	TH2F *fHistMultNonRemovedK0s; // enter comment here 
	TH2F *fHistMultNonRemovedLambdas; // enter comment here 
	TH2F *fHistMultNonRemovedElectrons; // enter comment here 
	TH2F *fHistMultNonRemovedPositrons; // enter comment here 
	TH2F *fHistMultNonRemovedMuPlus; // enter comment here
	TH2F *fHistMultNonRemovedMuMinus; // enter comment here
	TH2F *fHistMultNonRemovedNeutrons; // enter comment here
	TH2F *fHistMultNonRemovedAntiNeutrons; // enter comment here
	TH2F *fHistMultNonRemovedGammas; // enter comment here
	
	TH2F *fHistMultRemovedGammas; // enter comment here
	TH2F *fHistMultRemovedNeutrons; // enter comment here
	TH2F *fHistMultRemovedAntiNeutrons; // enter comment here
	
	TH2F *fHistTrackMultvsNonRemovedCharged; // enter comment here
	TH2F *fHistTrackMultvsNonRemovedNeutral; // enter comment here
	TH2F *fHistTrackMultvsRemovedGamma; // enter comment here
	
	TH2F *fHistClusterMultvsNonRemovedCharged; // enter comment here
	TH2F *fHistClusterMultvsNonRemovedNeutral; // enter comment here
	TH2F *fHistClusterMultvsRemovedGamma; // enter comment here
	
	TH2F *fHistMultvsNonRemovedChargedE; // enter comment here
	TH2F *fHistMultvsNonRemovedNeutralE; // enter comment here
	TH2F *fHistMultvsRemovedGammaE; // enter comment here
	
	Float_t fEtNonRemovedProtons; // enter comment here 
	Float_t fEtNonRemovedAntiProtons; // enter comment here 
	Float_t fEtNonRemovedPiPlus; // enter comment here 
	Float_t fEtNonRemovedPiMinus; // enter comment here 
	Float_t fEtNonRemovedKaonPlus; // enter comment here 
	Float_t fEtNonRemovedKaonMinus; // enter comment here 
	Float_t fEtNonRemovedK0s; // enter comment here 
	Float_t fEtNonRemovedLambdas; // enter comment here 
	Float_t fEtNonRemovedElectrons; // enter comment here 
	Float_t fEtNonRemovedPositrons; // enter comment here 
	Float_t fEtNonRemovedMuMinus; // enter comment here
	Float_t fEtNonRemovedMuPlus; // enter comment here
	Float_t fEtNonRemovedGammas; // enter comment here
	Float_t fEtNonRemovedGammasFromPi0; // enter comment here
	Float_t fEtNonRemovedNeutrons; // enter comment here
	Float_t fEtNonRemovedAntiNeutrons; // enter comment here
	
	Float_t fEtRemovedGammas; // enter comment here
	Float_t fEtRemovedNeutrons; // enter comment here
	Float_t fEtRemovedAntiNeutrons; // enter comment here
		
	Int_t fMultNonRemovedProtons; // enter comment here 
	Int_t fMultNonRemovedAntiProtons; // enter comment here 
	Int_t fMultNonRemovedPiPlus; // enter comment here 
	Int_t fMultNonRemovedPiMinus; // enter comment here 
	Int_t fMultNonRemovedKaonPlus; // enter comment here 
	Int_t fMultNonRemovedKaonMinus; // enter comment here 
	Int_t fMultNonRemovedK0s; // enter comment here 
	Int_t fMultNonRemovedLambdas; // enter comment here 
	Int_t fMultNonRemovedElectrons; // enter comment here 
	Int_t fMultNonRemovedPositrons; // enter comment here 
	Int_t fMultNonRemovedMuMinus; // enter comment here
	Int_t fMultNonRemovedMuPlus; // enter comment here
	Int_t fMultNonRemovedGammas; // enter comment here
	Int_t fMultNonRemovedNeutrons; // enter comment here
	Int_t fMultNonRemovedAntiNeutrons; // enter comment here
	
	Int_t fMultRemovedGammas; // enter comment here
	Int_t fMultRemovedNeutrons; // enter comment here
	Int_t fMultRemovedAntiNeutrons; // enter comment here
	
	Int_t fTrackMultInAcc; // enter comment here
	
	
	TH2F *fHistDxDzNonRemovedCharged; // enter comment here
	TH2F *fHistDxDzRemovedCharged; // enter comment here
	TH2F *fHistDxDzNonRemovedNeutral; // enter comment here
	TH2F *fHistDxDzRemovedNeutral; // enter comment here
	
	TH1F *fHistPiPlusMult; // enter comment here
	TH1F *fHistPiMinusMult; // enter comment here
	TH1F *fHistPiZeroMult; // enter comment here

	TH1F *fHistPiPlusMultAcc; // enter comment here
	TH1F *fHistPiMinusMultAcc; // enter comment here
	TH1F *fHistPiZeroMultAcc; // enter comment here
	
	Int_t fPiPlusMult; // enter comment here
	Int_t fPiMinusMult; // enter comment here
	Int_t fPiZeroMult; // enter comment here

	Int_t fPiPlusMultAcc; // enter comment here
	Int_t fPiMinusMultAcc; // enter comment here
	Int_t fPiZeroMultAcc; // enter comment here
	
	
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
