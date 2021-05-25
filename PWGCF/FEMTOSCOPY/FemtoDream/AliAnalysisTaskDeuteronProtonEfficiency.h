

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
    TH1F    *fHistPtProtonGen;
    TH1F    *fHistEtaProtonGen;
    TH1F    *fHistPtAntiProtonGen;
    TH1F    *fHistEtaAntiProtonGen;
    TH1F    *fHistPtDeuteronGen;
    TH1F    *fHistEtaDeuteronGen;
    TH1F    *fHistPtAntiDeuteronGen;
    TH1F    *fHistEtaAntiDeuteronGen;
    TH1F    *fHistPtHelium3Gen;
    TH1F    *fHistEtaHelium3Gen;
    TH1F    *fHistPtAntiHelium3Gen;
    TH1F    *fHistEtaAntiHelium3Gen;
    TH1F    *fHistSEDPairGen;
    TH1F    *fHistSEDAntiPairGen;
    TH2F    *fHistPtParticlesGen;
    TH2F    *fHistPtAntiParticlesGen;

    // Histograms for reconstructed particles
    TH1F    *fHistPtProtonRec;
    TH1F    *fHistEtaProtonRec;
    TH1F    *fHistPtAntiProtonRec;
    TH1F    *fHistEtaAntiProtonRec;
    TH1F    *fHistPtDeuteronRec;
    TH1F    *fHistEtaDeuteronRec;
    TH1F    *fHistPtAntiDeuteronRec;
    TH1F    *fHistEtaAntiDeuteronRec;
    TH1F    *fHistPtHelium3Rec;
    TH1F    *fHistEtaHelium3Rec;
    TH1F    *fHistPtAntiHelium3Rec;
    TH1F    *fHistEtaAntiHelium3Rec;
    TH1F    *fHistSEDPairRec;
    TH1F    *fHistSEDAntiPairRec;
    TH2F    *fHistPtParticlesRec;
    TH2F    *fHistPtAntiParticlesRec;

    TH1F    *fHistEventCounter;
    AliPIDResponse *fPIDResponse;

  ClassDef(AliAnalysisTaskDeuteronProtonEfficiency, 1); // analysisclass

}; // end of "public AliAnalysisTaskSE"

#endif
