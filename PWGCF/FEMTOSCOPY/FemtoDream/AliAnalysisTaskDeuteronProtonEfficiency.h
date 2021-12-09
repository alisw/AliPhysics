

#ifndef ALIANALYSISTASKDEUTERONPROTONEFFICIENCY_H 
#define ALIANALYSISTASKDEUTERONPROTONEFFICIENCY_H


#include "AliAnalysisTaskSE.h"



class AliESDEvent;
class AliESDInputHandler;


class AliAnalysisTaskDeuteronProtonEfficiency : public AliAnalysisTaskSE
{

  public:

    AliAnalysisTaskDeuteronProtonEfficiency();
    AliAnalysisTaskDeuteronProtonEfficiency(const char *name);
    AliAnalysisTaskDeuteronProtonEfficiency& operator = (const AliAnalysisTaskDeuteronProtonEfficiency &task);
    AliAnalysisTaskDeuteronProtonEfficiency(const AliAnalysisTaskDeuteronProtonEfficiency &task);
    virtual ~AliAnalysisTaskDeuteronProtonEfficiency();

    void UserCreateOutputObjects();
    void UserExec(Option_t *);
    void Terminate(Option_t *);
    double Bethe(const AliESDtrack &track, double mass, int charge, double *params);
    double CalculateRelativeMomentum(TLorentzVector &Pair, TLorentzVector &Part1, TLorentzVector &Part2);


  private:

    AliESDEvent		*fESD;
    AliESDInputHandler	*fESDHandler;
    AliVEvent		*fEvent;
    AliMCEvent		*mcEvent;
    bool		fMCtrue;
    AliESDtrackCuts	*fESDtrackCutsProton;
    AliESDtrackCuts	*fESDtrackCutsDeuteron;
    TList		*fHistList;

    // Histograms for generated particles

    // generated protons
    TH1F    *fHistPtProtonGenPDG;
    TH1F    *fHistPtProtonGenPrimary;
    TH1F    *fHistPtProtonGenEtaCut;
    TH1F    *fHistPtProtonGenPairPtCut;
    TH1F    *fHistEtaProtonGenPDG;
    TH1F    *fHistEtaProtonGenPrimary;
    TH1F    *fHistEtaProtonGenEtaCut;
    TH1F    *fHistEtaProtonGenPairPtCut;

    // generated antiprotons
    TH1F    *fHistPtAntiProtonGenPDG;
    TH1F    *fHistPtAntiProtonGenPrimary;
    TH1F    *fHistPtAntiProtonGenEtaCut;
    TH1F    *fHistPtAntiProtonGenPairPtCut;
    TH1F    *fHistEtaAntiProtonGenPDG;
    TH1F    *fHistEtaAntiProtonGenPrimary;
    TH1F    *fHistEtaAntiProtonGenEtaCut;
    TH1F    *fHistEtaAntiProtonGenPairPtCut;

    // generated deuterons
    TH1F    *fHistPtDeuteronGenPDG;
    TH1F    *fHistPtDeuteronGenPrimary;
    TH1F    *fHistPtDeuteronGenEtaCut;
    TH1F    *fHistPtDeuteronGenPairPtCut;
    TH1F    *fHistEtaDeuteronGenPDG;
    TH1F    *fHistEtaDeuteronGenPrimary;
    TH1F    *fHistEtaDeuteronGenEtaCut;
    TH1F    *fHistEtaDeuteronGenPairPtCut;

    // generated antideuterons
    TH1F    *fHistPtAntiDeuteronGenPDG;
    TH1F    *fHistPtAntiDeuteronGenPrimary;
    TH1F    *fHistPtAntiDeuteronGenEtaCut;
    TH1F    *fHistPtAntiDeuteronGenPairPtCut;
    TH1F    *fHistEtaAntiDeuteronGenPDG;
    TH1F    *fHistEtaAntiDeuteronGenPrimary;
    TH1F    *fHistEtaAntiDeuteronGenEtaCut;
    TH1F    *fHistEtaAntiDeuteronGenPairPtCut;

    // generated pairs
    TH1F    *fHistPtHelium3GenPairPtCut;
    TH1F    *fHistEtaHelium3GenPairPtCut;
    TH2F    *fHistPtParticlesGen;
    TH1F    *fHistSEDPairGen;

    // generated antipairs
    TH1F    *fHistPtAntiHelium3GenPairPtCut;
    TH1F    *fHistEtaAntiHelium3GenPairPtCut;
    TH2F    *fHistPtAntiParticlesGen;
    TH1F    *fHistSEDAntiPairGen;

    // reconstructed protons
    TH1F    *fHistPtProtonRecPDG;
    TH1F    *fHistPtProtonRecPrimary;
    TH1F    *fHistPtProtonRecTrackCuts;
    TH1F    *fHistPtProtonRecPairPtCut;
    TH1F    *fHistEtaProtonRecPDG;
    TH1F    *fHistEtaProtonRecPrimary;
    TH1F    *fHistEtaProtonRecTrackCuts;
    TH1F    *fHistEtaProtonRecPairPtCut;

    // reconstructed antiprotons
    TH1F    *fHistPtAntiProtonRecPDG;
    TH1F    *fHistPtAntiProtonRecPrimary;
    TH1F    *fHistPtAntiProtonRecTrackCuts;
    TH1F    *fHistPtAntiProtonRecPairPtCut;
    TH1F    *fHistEtaAntiProtonRecPDG;
    TH1F    *fHistEtaAntiProtonRecPrimary;
    TH1F    *fHistEtaAntiProtonRecTrackCuts;
    TH1F    *fHistEtaAntiProtonRecPairPtCut;

    // reconstructed deuterons
    TH1F    *fHistPtDeuteronRecPDG;
    TH1F    *fHistPtDeuteronRecPrimary;
    TH1F    *fHistPtDeuteronRecTrackCuts;
    TH1F    *fHistPtDeuteronRecPairPtCut;
    TH1F    *fHistEtaDeuteronRecPDG;
    TH1F    *fHistEtaDeuteronRecPrimary;
    TH1F    *fHistEtaDeuteronRecTrackCuts;
    TH1F    *fHistEtaDeuteronRecPairPtCut;

    // reconstructed antideuterons
    TH1F    *fHistPtAntiDeuteronRecPDG;
    TH1F    *fHistPtAntiDeuteronRecPrimary;
    TH1F    *fHistPtAntiDeuteronRecTrackCuts;
    TH1F    *fHistPtAntiDeuteronRecPairPtCut;
    TH1F    *fHistEtaAntiDeuteronRecPDG;
    TH1F    *fHistEtaAntiDeuteronRecPrimary;
    TH1F    *fHistEtaAntiDeuteronRecTrackCuts;
    TH1F    *fHistEtaAntiDeuteronRecPairPtCut;

    // reconstructed pairs
    TH1F    *fHistPtHelium3RecPairPtCut;
    TH1F    *fHistEtaHelium3RecPairPtCut;
    TH2F    *fHistPtParticlesRec;
    TH1F    *fHistSEDPairRec;

    // reconstructed antipairs
    TH1F    *fHistPtAntiHelium3RecPairPtCut;
    TH1F    *fHistEtaAntiHelium3RecPairPtCut;
    TH2F    *fHistPtAntiParticlesRec;
    TH1F    *fHistSEDAntiPairRec;

    // reconstructed particles and antiparticles
    TH2F    *fHistdEdx_LHC18a2a;
    TH2F    *fHistdEdx_LHC18a2b;
    TH2F    *fHistdEdx_LHC18a2b4;
    TH2F    *fHistdEdx_LHC20l7a;

    std::vector<int>*    GeneratedProtonArray;
    std::vector<int>*    GeneratedDeuteronArray;
    std::vector<int>*    GeneratedAntiProtonArray;
    std::vector<int>*    GeneratedAntiDeuteronArray;

    std::vector<int>*    ReconstructedProtonArray;
    std::vector<int>*    ReconstructedDeuteronArray;
    std::vector<int>*    ReconstructedAntiProtonArray;
    std::vector<int>*    ReconstructedAntiDeuteronArray;

    // event counter
    TH1F    *fHistEventCounter;
    AliPIDResponse *fPIDResponse;

  ClassDef(AliAnalysisTaskDeuteronProtonEfficiency, 1); // analysisclass

}; // end of "public AliAnalysisTaskSE"

#endif
