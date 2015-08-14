#ifndef ALIANALYSISTASKPIDconfig_H
#define ALIANALYSISTASKPIDconfig_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskPIDconfig.h */
// Author: Naghmeh Mohammadi, 10/07/2014

//==============================================================================
//
//
//
//
//==============================================================================

#include <TVectorDfwd.h>

#include "AliAnalysisTaskSE.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"
#include "TCutG.h"
#include "TF2.h"
#include "THnSparse.h"




class AliESDEvent;
class AliAODEvent;
class AliESDTrack;
class AliAODTrack;
class AliPIDResponse;
class TList;
class AliVEvent;
class TH1F;
class TH2F;
class TH3F;
class TH3F;



class AliAnalysisTaskPIDconfig : public AliAnalysisTaskSE {
    
    
public:
    AliAnalysisTaskPIDconfig();
    AliAnalysisTaskPIDconfig(const char *name);
    virtual ~AliAnalysisTaskPIDconfig();
    
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t * /*option*/);
    
    void SetFilterBit(Double_t b){fFilterBit = b;}
    void SetCentralityPercentileMin(Int_t b){fCentralityPercentileMin = b;}
    void SetCentralityPercentileMax(Int_t b){fCentralityPercentileMax = b;}
    void SetCentralityEstimator(TString b){fCentralityEstimator = b;}
    void SetUseCentrality(Bool_t b=kTRUE){fUseCentrality = b;}
    void SetCentralityTrigger(Int_t b=AliVEvent::kMB){ fTriggerSelection = b;}
    void SetDCAxyCut(Int_t b){fDCAxyCut = b;}
    void SetDCAzCut(Int_t b){fDCAzCut = b;}
    void SetCutTPCmultiplicityOutliersAOD(Bool_t b){fCutTPCmultiplicityOutliersAOD = b;}
    void SetData2011(Bool_t b){fData2011 = b;}
    void CheckCentrality(AliVEvent *event,Double_t centrality, Bool_t &centralitypass); //to use only events with the correct centrality....
    void SetCuts(Bool_t b){fPIDcuts = b;}
    //void MultiplicityOutlierCut(AliVEvent *event,Bool_t &centralitypass,Int_t ntracks);
    //void SetPIDcontoursList(TDirectory* b){fCutContourList = b;}
    
    void SetPIDPurityFunctions(Float_t PurityLevel);
    
protected:
    
    
    
private:
    AliVEvent             *fVevent;             //! event
    AliESDEvent           *fESD;                //! esd
    AliAODEvent           *fAOD;                //! aod
    AliPIDResponse        *fPIDResponse;        //! PID response Handler
    Int_t                  fTriggerSelection;   // trigger selection
    Int_t                  fCentralityPercentileMin;    // min centrality
    Int_t                  fCentralityPercentileMax;    // max
    Float_t                fPurityLevel;  //purity level
    Double_t               fFilterBit;                  // filterbit
    Double_t               fDCAxyCut;           // dca cut
    Double_t               fDCAzCut;            // dcz z
    Double_t               fLowPtPIDTPCnsigLow_Pion[6]; //nsig low cut Pion
    Double_t               fLowPtPIDTPCnsigHigh_Pion[6]; //nsig high cut Pion
    Double_t               fLowPtPIDTPCnsigLow_Kaon[6]; //nsig low cut Kaon
    Double_t               fLowPtPIDTPCnsigHigh_Kaon[6]; //nsig high cut Kaon
    Double_t               fvalueAll[4];//!All species info
    Bool_t                 fData2011;           // 2011 data
    Bool_t                 fUseCentrality;      // use centrality
    Bool_t                 fCutTPCmultiplicityOutliersAOD;      // do outlier cut
    Bool_t                 fPIDcuts;            // pid cuts
    TString                fCentralityEstimator;// cen estimator "V0M","TRK","TKL","ZDC","FMD"
    TFile                 *fPurityFunctionsFile;       //! purity functions file
    TDirectory            *fPurityFunctionsList;     //! purity functions list
    TList                 *fListQA;             //! List of all lists
    TList                 *fListQAtpctof;       //! List with combined PID from TPC + TOF
    TList                 *fListQAInfo;         //! list q ainfo
    TH1F                  *fhistCentralityPassBefore; //! cen histo before
    TH1F                  *fhistCentralityPassAfter; //! cen histo after
    TH1F                  *fNoEvents;           //! event no
    TH1F                  *fpVtxZ;              //! v vertex no
    TH2F                  *fhistDCABefore;      //! dca after hist
    TH2F                  *fhistDCAAfter;       //! another hist
    TH1F                  *fhistPhiDistBefore;  //! another hist
    TH1F                  *fhistPhiDistAfter;   //! another hist
    TH1F                  *fhistEtaDistBefore;  //! another hist
    TH1F                  *fhistEtaDistAfter;   //! another hist
    TH2F                  *fTPCvsGlobalMultBeforeOutliers;      //! another hist
    TH2F                  *fTPCvsGlobalMultAfterOutliers;       //! another hist
    TH2F                  *fTPCvsGlobalMultAfter;       //! another hist
    TH2F                  *fHistBetavsPTOFbeforePID;    //! another hist
    TH2F                  *fHistdEdxvsPTPCbeforePID;    //! another hist
    TH3F                  *fhistNsigmaP;        //! another hist
    TH2F                  *fhistTPCnSigmavsP;   //! another hist
    TH2F                  *fhistTOFnSigmavsP;   //! another hist
    TH2F                  *fHistBetavsPTOFafterPID;     //! another hist
    TH2F                  *fHistdEdxvsPTPCafterPID;     //! another hist
    TH2F                  *fHistBetavsPTOFafterPID_2;     //! another hist
    TH2F                  *fHistdEdxvsPTPCafterPID_2;     //! another hist
    TH2F                  *fHistBetavsPTOFafterPIDTPCTOF; //! another hist
    TH2F                  *fHistdEdxvsPTPCafterPIDTPCTOF; //! another hist
    TH2F                  *fHistBetavsPTOFafterPIDTPConly; //! another hist
    TH2F                  *fHistdEdxvsPTPCafterPIDTPConly; //! another hist
    TH2F                  *fHistPion_BetavsPTOFafterPIDTPCTOF; //! another hist
    TH2F                  *fHistPion_dEdxvsPTPCafterPIDTPCTOF; //!another hist
    TH2F                  *fHistKaon_BetavsPTOFafterPIDTPCTOF; //! another hist
    TH2F                  *fHistKaon_dEdxvsPTPCafterPIDTPCTOF; //!another hist
    TH2F                  *fHistProton_BetavsPTOFafterPIDTPCTOF; //! another hist
    TH2F                  *fHistProton_dEdxvsPTPCafterPIDTPCTOF; //!another hist
    TH1F                  *fhistPionEtaDistAfter; //!another hist
    TH1F                  *fhistKaonEtaDistAfter; //!another hist
    TH1F                  *fhistProtonEtaDistAfter; //!another hist
    THnSparseD            *fSparseAll;//!All species info
    TF2                   *fPurityFunction[180];  //purity functions
    
    //qa object initialisation
    void SetupTPCTOFqa();
    void SetupEventInfo();
    //
    
    AliAnalysisTaskPIDconfig(const AliAnalysisTaskPIDconfig &other);
    AliAnalysisTaskPIDconfig& operator=(const AliAnalysisTaskPIDconfig &other);
    
    ClassDef(AliAnalysisTaskPIDconfig,4)  // Task to properly set the PID response functions of all detectors
};

#endif
