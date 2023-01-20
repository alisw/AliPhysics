/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskStrangenessML_H
#define AliAnalysisTaskStrangenessML_H

#include "AliAnalysisTaskSE.h"
#include <TNamed.h>
#include "AliMachineLearning.h"
#include "AliNeuralNetwork.h"
#include "AliBDT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "TMath.h"
#include <iostream>
#include <TROOT.h>


class AliPIDResponse;

class AliAnalysisTaskStrangenessML : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskStrangenessML();
                                AliAnalysisTaskStrangenessML(const char *name);
        virtual                 ~AliAnalysisTaskStrangenessML();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

	//ML recipe
	void Set_NN_Recipe(TString FileName);
	void Set_BDT_Recipe(TString FileName);

	//Function to get cos(PA_BachelorBaryon)
	Float_t GetCosPA(AliESDtrack *lPosTrack, AliESDtrack *lNegTrack, AliESDEvent *lEvent);

    private:

        AliESDEvent*            fESD;           //! input event
        TList*                  fOutputList;    //! output list

	//For input ML Recipes
	AliNeuralNetwork* fXiMinusNN;
	AliBDT* fXiMinusBDT;
	
	//Output Histograms
	TH1F*			fHistCascadeCounts;
        TH1D*                   fHist_XiMinus_NNpred;
	TH3D*			fHist_XiMinus_NN_IM;
	TH1D*                   fHist_XiMinus_BDTpred;
        TH3D*                   fHist_XiMinus_BDT_IM;

	TH2F*			fHist_XiMinus_Std_IM;

	AliPIDResponse* 	fPIDResponse;	//! using other class here forward declared in [1]

	AliAnalysisTaskStrangenessML(const AliAnalysisTaskStrangenessML&); // not implemented
	AliAnalysisTaskStrangenessML& operator=(const AliAnalysisTaskStrangenessML&); // not implemented
        
	ClassDef(AliAnalysisTaskStrangenessML, 1);
};

#endif
