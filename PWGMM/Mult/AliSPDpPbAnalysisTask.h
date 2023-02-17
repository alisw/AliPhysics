#ifndef AliSPDpPbAnalysisTask_H
#define AliSPDpPbAnalysisTask_H

//pA data analysis (in construction)
//Author: Eliana Marroquin (eliana.marroquin@cern.ch)
//Inspired by Abhi Modak (abhi.modak@cern.ch) code for pPb collisions 

class TH1F;
class TH2F;
class TString;
class TNtuple;
class TArrayF;
class AliESDEvent;
class AliHeader;
class TParticle;
class AliMultSelection;
class AliVMultiplicity;
class AliMCEvent;
class AliVEvent;
class AliStack;
#include "AliAnalysisTaskSE.h"

class AliSPDpPbAnalysisTask : public AliAnalysisTaskSE
{
public:
    // two class constructors
    AliSPDpPbAnalysisTask();
    AliSPDpPbAnalysisTask(const char *name);

    // class destructor
    virtual ~AliSPDpPbAnalysisTask();

    // Methods of the analysis framework
    // called once at the beginning of the runtime
    virtual void UserCreateOutputObjects();
    // called for each event
    virtual void UserExec(Option_t *);
    // called at end of analysis
    virtual void Terminate(Option_t *);
    void SetTrigger(const char* trigger){
    fTrigSel = trigger;
    }
    void SetCentralityEstimator(const char* Estimator){
    fCentEstimator = Estimator;
    }

private:

    AliESDEvent                 *fESD;                         //! input event
    AliMCEvent                  *fMCEvent;                     //! corresponding MC event
    AliStack                    *fMCStack;
    AliESDVZERO                 *fEsdV0;
    TList                       *fOutputList;                  //! output list
    TH1F                        *fHistVtxZ;                    //! primary vertex Z 
    TH1F                        *fHistTotEvent;                //! total event
    TH1F                        *fNTracklets;                  //! number of tracklets per event
    TH2F                        *fHistTrksVtxZ;
    TH2F                        *fHistVtxXY;
    TH2F                        *fHistVtxZEta;
    TH2F                        *fHistEtaVtxZ;
    TH2F                        *fHistPhiEta;
    TH2F                        *fHistEtaPhi;
    TH1F                        *fHistEta;
    TH1F                        *fHistMult;
    TH1F                        *fHistTrks05;
    TH1F                        *fHistTrks10;
    TH1F                        *fHistTrks15;
    TH1F                        *fMCHistVtxZ;
    TH2F                        *fMCHistXY;
    TH2F                        *fMCHistPhiEta;
    TH2F                        *fMCHistPhiEta10;
    TH2F                        *fMCHistEtaVxtZ;
    TH1F                        *fMCHistPt;
    TH1F                        *fMCHistEta;
    TH1F                        *fMCHistEta10;
    TH1F                        *fMCHistPhi;
    TH1F                        *fMCHistNch;
    TH1F                        *fMCHistNch10;
    TH1F                        *fMCHistNch05;
    TH1F                        *fMCHistNch15;
    TH2F                        *fResponseMatrix;
    TH2F                        *fResponseMatrix10;
    TH2F                        *fResponseMatrix101;
    TH2F                        *fResponseMatrix05;
    TH2F                        *fResponseMatrix051;
    TH2F                        *fResponseMatrix15;
    TH2F                        *fResponseMatrix151;
    TH2F                        *fEffVsNch;
    TString                     fCentEstimator;
    TString                     fTrigSel;
    AliMultSelection*           MultSelection=nullptr;
    AliVMultiplicity*           fMultiplicity=nullptr;


    AliSPDpPbAnalysisTask(const AliSPDpPbAnalysisTask&); // not implemented
    AliSPDpPbAnalysisTask& operator=(const AliSPDpPbAnalysisTask&); // not implemented
    
    ClassDef(AliSPDpPbAnalysisTask, 1);
};

#endif