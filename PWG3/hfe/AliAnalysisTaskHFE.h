/**************************************************************************
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
#ifndef ALIANALYSISTASKHFE_H
#define ALIANALYSISTASKHFE_H

#include "AliAnalysisTask.h"

class AliHFEpid;
class AliHFEcuts;
class AliCFManager;
class AliESDEvent;
class AliESDtrackCuts;
class AliMCEvent;
class TH1I; 
class TList;

class AliAnalysisTaskHFE : public AliAnalysisTask{
  enum{
    kIsQAOn = BIT(14)
  };
	public:
		AliAnalysisTaskHFE();
		~AliAnalysisTaskHFE();

		virtual void ConnectInputData(Option_t *);
		virtual void CreateOutputObjects();
		virtual void Exec(Option_t *);
		virtual void Terminate(Option_t *);

    void SetQAOn() { SetBit(kIsQAOn, kTRUE); };
    Bool_t IsQAOn() const { return TestBit(kIsQAOn); };

	private:
    void MakeParticleContainer();

		AliESDEvent *fESD;		                //! The ESD Event
		AliMCEvent *fMC;		                  //! The MC Event
		AliCFManager *fCFM;		                //! Correction Framework Manager
		AliHFEpid *fPID;         //! PID
    AliHFEcuts *fCuts;     //! Cut Collection
		TH1I *fNEvents;			                  //! counter for the number of Events
		TList *fQA;			                      //! QA histos for the cuts

	ClassDef(AliAnalysisTaskHFE, 1)    // The electron Analysis Task
};
#endif
