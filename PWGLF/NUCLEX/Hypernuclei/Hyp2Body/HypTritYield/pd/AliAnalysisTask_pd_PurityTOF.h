#ifndef AliAnalysisTask_pd_PurityTOF_H
#define AliAnalysisTask_pd_PurityTOF_H

#include "AliAnalysisTaskSE.h"
#include "TObject.h"
#include "AliEventCuts.h"

class AliAODEvent;
class AliAODInputHandler;







class AliAnalysisTask_pd_PurityTOF : public AliAnalysisTaskSE
{

  public:

    AliAnalysisTask_pd_PurityTOF();
    AliAnalysisTask_pd_PurityTOF(const char *name,Int_t CollisionSystem, Bool_t UseOpenCuts, Bool_t isMC, Bool_t SavePairsOnly);
    AliAnalysisTask_pd_PurityTOF& operator = (const AliAnalysisTask_pd_PurityTOF &task);
    AliAnalysisTask_pd_PurityTOF(const AliAnalysisTask_pd_PurityTOF &task);
    virtual ~AliAnalysisTask_pd_PurityTOF();

    void UserCreateOutputObjects();
    void UserExec(Option_t *);
    void Terminate(Option_t *);
    Double_t CalculateBetaTOF(AliAODTrack &track); 
    Double_t CalculateMassSquareTOF(AliAODTrack &track);
    Double_t CalculateSigmaMassSquareTOF(Double_t pT, Double_t massSq, Int_t ParticleSpecies, Int_t RunNumber);
    Bool_t CheckProtonCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, Int_t ParticleSpecies, Int_t RunNumber);
    Bool_t CheckDeuteronCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, Int_t ParticleSpecies, Int_t RunNumber);
    Double_t CalculateSigmadEdxTPC(AliAODTrack &Track, Int_t ParticleSpecies);
    Double_t CalculateSigmadEdxITS(AliAODTrack &Track, Int_t ParticleSpecies, Int_t RunNumber);


  private:

    AliAODEvent		  *fAODEvent;
    AliAODInputHandler	  *fAODHandler; 
    AliAODHeader	  *fHeader;
    AliPIDResponse	  *fPIDResponse;

    Int_t   fCollisionSystem;
    Bool_t  fUseOpenCuts;
    Bool_t  fIsMC;
    Bool_t  fSavePairsOnly;
    AliTimeRangeCut fTimeRangeCut;

    TList *OutputList;
    TH2F *h_Purity_Proton;
    TH2F *h_Purity_AntiProton;
    TH2F *h_Purity_Deuteron;
    TH2F *h_Purity_AntiDeuteron;
    TH2F *h_Purity2_Deuteron;
    TH2F *h_Purity2_AntiDeuteron;



    ClassDef(AliAnalysisTask_pd_PurityTOF,1);

};






#endif
