/*************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//

// Author:
// Paul Buehler <paul.buehler@oeaw.ac.at>

#ifndef ALIANALYSISTASKCEPQA_H
#define ALIANALYSISTASKCEPQA_H

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskCEPQA : public AliAnalysisTaskSE
{
	
	public:

		AliAnalysisTaskCEPQA(const char* name);
		AliAnalysisTaskCEPQA();
		virtual ~AliAnalysisTaskCEPQA();

		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t *);
		virtual void Terminate(Option_t *);
		virtual void Clear(Option_t *);

	private:

		AliAnalysisTaskCEPQA(const AliAnalysisTaskCEPQA &p);
		AliAnalysisTaskCEPQA& operator=(const AliAnalysisTaskCEPQA &p);

		// functions used in UserExec
		void CheckVertex    ();
		void CheckSPDactive ();
		void CheckV0active  ();
		void CheckFMDactive ();
		void CheckADactive  ();
		void CheckZDNactive ();
		void CheckZDPactive ();
		void CheckGapType   ();
		Bool_t FilterEvent  ();

		Bool_t CheckInput();
		void PostOutputs();
		
		// Variables
    TFile *ff;
    AliESDEvent *fESDEvent;
    AliTriggerAnalysis* fTrigger;

		// Output objects
		TTree *fTree;
		TList *fList;
    
		TH1D *fHistStatus;

		Int_t    fRunNum;
		Int_t   fEvCounter;
		Bool_t   fisMBOR;
		Bool_t   fisMBAND;
		Bool_t   fisPileup;
		Int_t    fSPDmul;
		Int_t    fVtype;
		Double_t fVposx;
		Double_t fVposy;
		Double_t fVposz;
		Bool_t   fisV0A;
		Bool_t   fisV0C;
		Int_t    fFMDAnumHC;
		Int_t    fFMDCnumHC;
		Float_t  fADAmul;
		Float_t  fADCmul;
		Int_t    fADATCh;
		Int_t    fADCTCh;
		Bool_t   fADATD;
		Bool_t   fADCTD;
		Float_t  fZDCTime;
		Double_t fZNAEne;
		Double_t fZNCEne;
		Double_t fZPAEne;
		Double_t fZPCEne;
    
    Bool_t   fisNG;
    Bool_t   fisSG;
    Bool_t   fisDG;
    		
		ClassDef(AliAnalysisTaskCEPQA, 1);

};

#endif
