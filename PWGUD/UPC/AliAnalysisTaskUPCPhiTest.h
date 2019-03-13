/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskUPCPhiTest_H
#define AliAnalysisTaskUPCPhiTest_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskUPCPhiTest : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskUPCPhiTest();
                                AliAnalysisTaskUPCPhiTest(const char *name);
        virtual                 ~AliAnalysisTaskUPCPhiTest();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        AliAODEvent*            fAOD;  //! input event
        
        AliPIDResponse*         fPIDResponse;            
        TList*                  fOutputList;    //! output list
        TList*                  fOutputList2;
        TTree*                  fTreeP_TPC;
        TTree*                  fTree_NoCut;
        //One Dimensional Counters
        TH1D* fHistRunCounter;    
        TH1I* fHistCounter;     
        TH1D* fHistCcup4Triggers;
        TH1D* fHistCcup7Triggers;
        TH1D* fHistCcup2Triggers;
        TH1D* fHistCint1Triggers;
        TH1D* fHistCint6Triggers;
        TH1D* fHistC0tvxAndCint1;
        TH1D* fHistZedTriggers;
        TH1D* fHistCvlnTriggers;
        TH1D* fHistMBTriggers;
        TH1D* fHistCentralTriggers;
        TH1D* fHistSemiCentralTriggers;
        
        TH1D* fHistCTest58Triggers;
        TH1D* fHistCTest59Triggers;
        TH1D* fHistCTest60Triggers;
        TH1D* fHistCTest61Triggers;
        
        TH1D* fHistCcup8Triggers;
        TH1D* fHistCcup9Triggers;
        TH1D* fHistCcup291Triggers;
        TH1D* fHistCcup301Triggers;
        TH1D* fHistCcup311Triggers;
        TH1D* fHistCcup29Triggers;
        TH1D* fHistCcup30Triggers;
        TH1D* fHistCcup31Triggers;
        TH1D* fHistCtrueTriggers;
        // information that are being saved
        Float_t                 fPt;
        Float_t                 fPt0; //pt of daughter 1
        Float_t                 fPt1; // pt of daughter 2
        Int_t                   fRunNumber; 
        Float_t                 fPhi1;
        Float_t                 fPhi2;
        Float_t                 fDelPhi;
        Int_t                  fITSInHits;
        Int_t                  fITSOutHits;
        
        //parent mass calculated considering differnt daughter mass
        
        Float_t                 fM; 
        Float_t                 fPiM;
        Float_t                 fMuM;
        Float_t                 fElM;
       //sigma signal from TPC
        Float_t                 fKaonSigma0;
        Float_t                 fKaonSigma1;
        Float_t                 fMuSigma0;
        Float_t                 fMuSigma1;
        Float_t                 fElSigma0;
        Float_t                 fElSigma1;
        Float_t                 fPiSigma0;
        Float_t                 fPiSigma1;
        //Sigma signal from TOF
        Float_t                 fKaonSigmaTOF0;
        Float_t                 fKaonSigmaTOF1;
        Float_t                 fPiSigmaTOF0;
        Float_t                 fPiSigmaTOF1;
        Float_t                 fMuSigmaTOF0;
        Float_t                 fMuSigmaTOF1;
        Float_t                 fElSigmaTOF0;
        Float_t                 fElSigmaTOF1;
        
        //differnt informations  1 is for daughter 1 and 2 is for daughte 2
        Float_t                 fTPCcluster1;
        Float_t                 fEta1;
        Float_t                 fTPCcluster2;
        Float_t                 fEta2;
        Float_t                 fDCAxy1;
        Float_t                 fDCAz1;
        Float_t                 fDCAxy2;
        Float_t                 fDCAz2;
        //saving trigger inforamation so that different trigger condition can be observed individucally
        Int_t                   fTriggerClass;
        Float_t                 fPd0;
        Float_t                 fPd1;
        Float_t                 fPp;
        
        Float_t                 fdEdX0;
        Float_t                 fdEdX1;
        
        Float_t                 fdEdXTOF0;
        Float_t                 fdEdXTOF1; 
        //Storing Charge
	Int_t			fCharge0;
	Int_t			fCharge1;
        Float_t                 fZNAenergy; 
        Float_t                 fZNCenergy;     
         
        Float_t                 fZDCAtime;
        Float_t                 fZDCCtime;
        
        AliAnalysisTaskUPCPhiTest(const AliAnalysisTaskUPCPhiTest&); // not implemented
        AliAnalysisTaskUPCPhiTest& operator=(const AliAnalysisTaskUPCPhiTest&); // not implemented

        ClassDef(AliAnalysisTaskUPCPhiTest, 1);
};

#endif
