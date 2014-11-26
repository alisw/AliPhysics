/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*
 * A small task dumping all EMCal trigger related information into a TTree
 *      Author: Markus Fasel
 */

#ifndef ALIANALYSISTASKEMCALTRIGGERTREEWRITER_H_
#define ALIANALYSISTASKEMCALTRIGGERTREEWRITER_H_

#include "AliAnalysisTaskSE.h"

class TTree;

class AliAnalysisTaskEmcalTriggerTreeWriter : public AliAnalysisTaskSE {
public:
	AliAnalysisTaskEmcalTriggerTreeWriter();
	AliAnalysisTaskEmcalTriggerTreeWriter(const char *name);
	virtual ~AliAnalysisTaskEmcalTriggerTreeWriter();

	virtual void UserCreateOutputObjects();
	virtual void UserExec(Option_t *);

private:
	struct TriggerInfo{
		Int_t fRun;
		Int_t fCol;
		Int_t fRow;
		Int_t fNL0Times;
		Int_t fLevel0Times[10];
		Int_t fADC;
		Float_t fAmplitude;
		Float_t fTime;
		Int_t fTriggerBits;
		Int_t fL1Threshold;
		Int_t fL1V0;

		TriggerInfo():
			fRun(0),
			fCol(0), fRow(0),
			fNL0Times(0), fADC(0), fAmplitude(0.),
			fTime(0), fTriggerBits(0), fL1Threshold(0), fL1V0(0)
		{
			memset(fLevel0Times, 0, sizeof(Int_t) * 10);
		}
		void Reset(){
			fRun = 0; fCol = 0; fRow = 0;
			fNL0Times = 0; fADC = 0; fAmplitude = 0.;
			fTime = 0.; fTriggerBits = 0; fL1Threshold = 0; fL1V0 = 0;
			memset(fLevel0Times, 0, sizeof(Int_t) * 10);
		}
	};
	TTree *fOutputTree;						//! Output tree with tracks
	TriggerInfo fOutputInfo;					// Track Info for the tree

	ClassDef(AliAnalysisTaskEmcalTriggerTreeWriter, 1)
};

#endif /* ALIANALYSISTASKEMCALTRIGGERTREEWRITER_H_ */
