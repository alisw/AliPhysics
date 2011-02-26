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
	
 private:

    //Declare it private to avoid compilation warning
    AliAnalysisEtMonteCarlo & operator = (const AliAnalysisEtMonteCarlo & g) ;//cpy assignment
    AliAnalysisEtMonteCarlo(const AliAnalysisEtMonteCarlo & g) ; // cpy ctor
    ClassDef(AliAnalysisEtMonteCarlo, 2);
};

#endif // ALIANALYSISETMONTECARLO_H
