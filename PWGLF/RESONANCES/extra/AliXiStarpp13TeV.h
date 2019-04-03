#ifndef ALIXISTARPP13TEV_H
#define ALIXISTARPP13TEV_H
//
// Class AliXiStarpp13TeV
//
// AliXiStarpp13TeV
// author:
//  (Original Code) Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//  (1st Modification) Jihye Song (jihye.song@cern.ch)
//  (last Modification) Bong-Hwi Lim (bong-hwi.lim@cern.ch)


class TH1F;
class TH1D;
class TH2D;
class TH3D;
class TProfile;
class TTree;

class AliESDEvent;
class AliESDtrackCuts;
class AliESDpid;

#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODPid.h"
#include "AliESDpid.h"
#include "AliXiStarpp13TeVEventCollection.h"
#include "AliESDVZERO.h"
#include "AliESDTZERO.h"
#include "AliVertex.h"


class AliXiStarpp13TeV : public AliAnalysisTaskSE {
public:
    
    AliXiStarpp13TeV();
    AliXiStarpp13TeV(const char *name, Bool_t AODdecision, Bool_t MCdecision, Int_t CutListOption=0, Bool_t DevelopmentMode = kFALSE, Bool_t fHMTrigger = kFALSE, Bool_t fPIDOption = kFALSE, Bool_t SetSystematic = kTRUE);

    virtual ~AliXiStarpp13TeV();
    AliXiStarpp13TeV(const AliXiStarpp13TeV &obj );
    AliXiStarpp13TeV &operator=(const AliXiStarpp13TeV &obj );

    enum {
        kNbinsM              = 200, // mult bins for certain histograms //300
        kXiStarCode          = 3324,// Xi(1530)^0 MC code
        kXiCode              = 3312,// Xi- MC code
        kLambdaCode          = 3122,// Lambda MC code
        kProtonCode          = 2212,// Proton+ MC code
        kPionCode            = 211,// Pion+ MC code
        kNCutVariations      = 21,// number of cut variations // 13
        kNCuts               = 13// number of cut types //15
    };

    //=================================================================================//
    //generated Histograms//
    //=================================================================================//
    
    struct St_CutType {
        TH3F *fXi; //!
        TH3F *fXibar; //!
        //
        TH3F *fXiMinusPiPlus; //!
        TH3F *fXiMinusPiMinus; //!
        TH3F *fXiPlusPiPlus; //!
        TH3F *fXiPlusPiMinus; //!
        
        TH3F *fXiMinusPiPlusbkg; //!
        TH3F *fXiMinusPiMinusbkg; //!
        TH3F *fXiPlusPiPlusbkg; //!
        TH3F *fXiPlusPiMinusbkg; //!
        //
        TH3F *fMCrecXi; //!
        TH3F *fMCrecXibar; //!
        
        TH3F *fMCrecXiMinusPiPlus; //!
        TH3F *fMCrecXiPlusPiMinus; //!
       
    };
    struct St_CutType CutVar[kNCutVariations]; //!

private:
    
    virtual void   UserCreateOutputObjects();
    virtual void   Exec(Option_t *option);
    virtual void   Terminate(Option_t *);
    void XiStarInit();// initialization of fixed values
    Double_t LinearPropagateToDCA(AliESDtrack*, AliESDtrack*, Double_t);// for linear propagation
    Double_t Det(Double_t, Double_t, Double_t, Double_t) const;// for linear propagation
    Double_t Det(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t) const;// for linear propagation
     ULong64_t GetMCEventNumber();

    /*
    TH3F *fVertexDist1;
    TH3F *fVertexDist3;
    TH2F *fDCADist;
    TH3F *fMultDist3d;
    TH1F *fMultDist1;
    TH1F *fMultDist2;
    TH1F *fMultDist3;
    TH1F *fMultDist4;
    TH1F *fMultDist5;
    TH1F *fMultDist_pp;
    TH1F *fMultDist_pp_afterPileUpReject;
    TH1F *hEventSelecInfo;
    TH1F *hNumberOfEvent;

    TH1F *fPtDist;
    TH1F *fPhiDist;
    TH1F *fEtaDist;
    TH1F *fXiStarYDist;
    TH1F *fQAXiStarYDist;

    TH1F *fXiStarYDistMC;
    TH1F *fXiYDistMC1;
    TH1F *fXiStarYDistMC1;
    TH1F *fXiYDistMCout;
    TH1F *fXiStarYDistMCout;
    TH1F *fCutEvents;
    */

    
    const char* fname;// name of class
   // AliInputEventHandler *fEventHandler;                              //  for ESDs or AODs
    AliESDEvent            *fESD; //!    // ESD object
    TList                  *fOutputList; //! Compact Output list
    AliESDtrackCuts        *fTrackCut; //! ESD track cuts
    AliPIDResponse         *fPIDResponse; //! PID object
    
    AliXiStarpp13TeVEventCollection ***fEC; //! The event collection
    AliXiStarpp13TeVEventStruct *fEvt; //! The current event type
    AliXiStarpp13TeVTrackStruct *fTempStruct; //! A temporary track storage.  Eventually put into fEvt
    
    //
    
    Int_t fZvertexBins;// number of Z-vertex bins for event-mixing
    Int_t fEventsToMix;// number of maximum events to mix
    Int_t fMultBins;// number of multiplicity bins for event-mixing
    Int_t fMultLimits[11+1];// the multiplicity edges of the mult bins
    Bool_t fMCcase;// switch for MC data or real data
    Bool_t fAODcase;// switch for AODs or ESDs
    Bool_t fDevelopeMode;
    Bool_t fHMTrigger;
    Bool_t fPIDOption;
    Bool_t fSetSystematic;
    Int_t fEventCounter;// The event counter
    ULong64_t fEventNumber; // calcuate event number

    // cut list data members
    Float_t fMaxDecayLength;// max decay length
    Float_t fMassWindow;// Mass window of acceptance for Lambda and Xi candidates
    
    Double_t fCovMatrix[21];// Covarience matrix of track
    Double_t fTrueMassPr, fTrueMassPi, fTrueMassK, fTrueMassLam, fTrueMassXi;// The PDG mass values
    Bool_t IsTPC  (AliESDtrack *track);

  
    

    
    
    AliESDtrack* fESDTrack4; //! esdtrack for XiStar's daughter pion
    AliESDtrack* fXiTrack; //! esdtrack for XiStar's daughter Xi
    
    Int_t fCutList;// Cut List option (mean values or systematic variations)
    
    Float_t fDecayParameters[kNCuts];// array of reconstruction kinematics
    Float_t fCutValues[kNCutVariations][kNCuts];// array of reconstruction kinematics
    
    ClassDef(AliXiStarpp13TeV, 2); 
};

#endif