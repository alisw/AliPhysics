/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskTransTask_H
#define AliAnalysisTaskTransTask_H

#include "AliAnalysisTaskSE.h"

class TTree;
class TH1I;
class TList;
class AliAODEvent;

class AliAnalysisTaskTransTask : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskTransTask();
                                AliAnalysisTaskTransTask(const char *name);
        virtual                 ~AliAnalysisTaskTransTask();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
	TTree *fAnaTree; //! analysis tree
	Int_t fRunNum;
	Int_t fTracklets;
	Int_t fCtrue;
	Int_t fC1zed;	

	UInt_t fL0inputs;
	UInt_t fL1inputs;	

	Double_t fZem1Energy;
	Double_t fZem2Energy; 	

	Double_t fZNCEnergy; 
	Double_t fZNAEnergy;
	Double_t fZPCEnergy; 
	Double_t fZPAEnergy;
	Double_t fZNATDC[4];
	Double_t fZNCTDC[4];
	Double_t fZPATDC[4];
	Double_t fZPCTDC[4];
	Int_t fV0ADecision; 
	Int_t fV0CDecision;
	Int_t fADADecision; 
	Int_t fADCDecision;
  	TBits fIR1Map;
  	TBits fIR2Map;
	TH1I *fCounter; //! analysis counter
	
        AliAnalysisTaskTransTask(const AliAnalysisTaskTransTask&); // not implemented
        AliAnalysisTaskTransTask& operator=(const AliAnalysisTaskTransTask&); // not implemented

        ClassDef(AliAnalysisTaskTransTask, 1);
};

#endif
