#ifndef AliSPDpPbAnalysisTaskData_H
#define AliSPDpPbAnalysisTaskData_H

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
class AliVEvent;
#include "AliAnalysisTaskSE.h"

class AliSPDpPbAnalysisTaskData : public AliAnalysisTaskSE
{
public:
    // two class constructors
    AliSPDpPbAnalysisTaskData();
    AliSPDpPbAnalysisTaskData(const char *name);

    // class destructor
    virtual ~AliSPDpPbAnalysisTaskData();

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

    TString                     fCentEstimator;
    TString                     fTrigSel;
    AliMultSelection*           MultSelection=nullptr;
    AliVMultiplicity*           fMultiplicity=nullptr;


    AliSPDpPbAnalysisTaskData(const AliSPDpPbAnalysisTaskData&); // not implemented
    AliSPDpPbAnalysisTaskData& operator=(const AliSPDpPbAnalysisTaskData&); // not implemented
    
    ClassDef(AliSPDpPbAnalysisTaskData, 1);
};

#endif