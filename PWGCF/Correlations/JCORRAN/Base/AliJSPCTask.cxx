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

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for providing various flow informations
// author: D.J. Kim(dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla
// Finland
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////
#include <AliAnalysisManager.h>
#include <AliJBaseTrack.h>
#include "AliJSPCTask.h"

//______________________________________________________________________________
AliJSPCTask::AliJSPCTask() :   
  AliAnalysisTaskSE("JSPCTask"),
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fIsMC(kFALSE),
	fSPC(NULL),
	bJSPCSaveAllQA(kTRUE),
	fJSPCcent_0(0.), fJSPCcent_1(0.), fJSPCcent_2(0.), fJSPCcent_3(0.), fJSPCcent_4(0.), fJSPCcent_5(0.), fJSPCcent_6(0.), fJSPCcent_7(0.), fJSPCcent_8(0.), fJSPCcent_9(0.), 
  fJSPCcent_10(0.), fJSPCcent_11(0.), fJSPCcent_12(0.), fJSPCcent_13(0.), fJSPCcent_14(0.), fJSPCcent_15(0.), fJSPCcent_16(0.),
  bJSPCDoFisherYates(kFALSE),
  fJSPCFisherYatesCutOff(1.),
  fJSPCMinNumberPart(14),
  bJSPCUseWeightsNUE(kTRUE),
  bJSPCUseWeightsNUA(kFALSE),
  fJSPCNumber(0),  		//number of correlation first correlator
  fJSPCNumberSecond(0), 	//number of correlation second correlator
  fJSPCNumberThird(0),	//number of correlation second correlator
  fJSPCNumberFourth(0),	//number of correlation fourth correlator
  fJSPCNumberFifth(0),	//number of correlation fifth correlator
  fJSPCNumberSixth(0),	//number of correlation sixth correlator
  fJSPCNumberSeventh(0),	//number of correlation seventh correlator
  fJSPCNumberEighth(0),	//number of correlation eighth correlator
  fJSPCa1(0), fJSPCa2(0), fJSPCa3(0), fJSPCa4(0), fJSPCa5(0), fJSPCa6(0), fJSPCa7(0), //first set of harmonics
  fJSPCb1(0), fJSPCb2(0), fJSPCb3(0), fJSPCb4(0), fJSPCb5(0), fJSPCb6(0), fJSPCb7(0), //second set of harmonics
  fJSPCd1(0), fJSPCd2(0), fJSPCd3(0), fJSPCd4(0), fJSPCd5(0), fJSPCd6(0), fJSPCd7(0), //third set of harmonics
  fJSPCe1(0), fJSPCe2(0), fJSPCe3(0), fJSPCe4(0), fJSPCe5(0), fJSPCe6(0), fJSPCe7(0), //fourth set of harmonics
  fJSPCf1(0), fJSPCf2(0), fJSPCf3(0), fJSPCf4(0), fJSPCf5(0), fJSPCf6(0), fJSPCf7(0),  //fifth set of harmonics 
  fJSPCg1(0), fJSPCg2(0), fJSPCg3(0), fJSPCg4(0), fJSPCg5(0), fJSPCg6(0), fJSPCg7(0),  // sixth set of harmonics 
  fJSPCh1(0), fJSPCh2(0), fJSPCh3(0), fJSPCh4(0), fJSPCh5(0), fJSPCh6(0), fJSPCh7(0), //seventh set of harmonics
  fJSPCi1(0), fJSPCi2(0), fJSPCi3(0), fJSPCi4(0), fJSPCi5(0), fJSPCi6(0), fJSPCi7(0),  //eighth set of harmonics 
  bJSPCDoMixed(kFALSE),
  bJSPCDifferentCharge(kTRUE),
  bJSPCSetSameChargePositive(kTRUE),
  fJSPCMixedHarmonic(0), 
  bJSPCComputeEtaGap(kFALSE),		// Do eta gap computation if kTRUE. Default kFALSE
  fJSPCEtaGap(0.8)			// Value of eta gap
{
}

//______________________________________________________________________________
AliJSPCTask::AliJSPCTask(const char *name):
	AliAnalysisTaskSE(name), 
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fIsMC(kFALSE),
	fSPC(NULL),
	bJSPCSaveAllQA(kTRUE),
	fJSPCcent_0(0.), fJSPCcent_1(0.), fJSPCcent_2(0.), fJSPCcent_3(0.), fJSPCcent_4(0.), fJSPCcent_5(0.), fJSPCcent_6(0.), fJSPCcent_7(0.), fJSPCcent_8(0.), fJSPCcent_9(0.), 
  fJSPCcent_10(0.), fJSPCcent_11(0.), fJSPCcent_12(0.), fJSPCcent_13(0.), fJSPCcent_14(0.), fJSPCcent_15(0.), fJSPCcent_16(0.),
  bJSPCDoFisherYates(kFALSE),
  fJSPCFisherYatesCutOff(1.),
  fJSPCMinNumberPart(14),
  bJSPCUseWeightsNUE(kTRUE),
  bJSPCUseWeightsNUA(kFALSE),
  fJSPCNumber(0),  		//number of correlation first correlator
  fJSPCNumberSecond(0), 	//number of correlation second correlator
  fJSPCNumberThird(0),	//number of correlation second correlator
  fJSPCNumberFourth(0),	//number of correlation fourth correlator
  fJSPCNumberFifth(0),	//number of correlation fifth correlator
  fJSPCNumberSixth(0),	//number of correlation sixth correlator
  fJSPCNumberSeventh(0),	//number of correlation seventh correlator
  fJSPCNumberEighth(0),	//number of correlation eighth correlator
  fJSPCa1(0), fJSPCa2(0), fJSPCa3(0), fJSPCa4(0), fJSPCa5(0), fJSPCa6(0), fJSPCa7(0), //first set of harmonics
  fJSPCb1(0), fJSPCb2(0), fJSPCb3(0), fJSPCb4(0), fJSPCb5(0), fJSPCb6(0), fJSPCb7(0), //second set of harmonics
  fJSPCd1(0), fJSPCd2(0), fJSPCd3(0), fJSPCd4(0), fJSPCd5(0), fJSPCd6(0), fJSPCd7(0), //third set of harmonics
  fJSPCe1(0), fJSPCe2(0), fJSPCe3(0), fJSPCe4(0), fJSPCe5(0), fJSPCe6(0), fJSPCe7(0), //fourth set of harmonics
  fJSPCf1(0), fJSPCf2(0), fJSPCf3(0), fJSPCf4(0), fJSPCf5(0), fJSPCf6(0), fJSPCf7(0),  //fifth set of harmonics 
  fJSPCg1(0), fJSPCg2(0), fJSPCg3(0), fJSPCg4(0), fJSPCg5(0), fJSPCg6(0), fJSPCg7(0),  // sixth set of harmonics 
  fJSPCh1(0), fJSPCh2(0), fJSPCh3(0), fJSPCh4(0), fJSPCh5(0), fJSPCh6(0), fJSPCh7(0), //seventh set of harmonics
  fJSPCi1(0), fJSPCi2(0), fJSPCi3(0), fJSPCi4(0), fJSPCi5(0), fJSPCi6(0), fJSPCi7(0),  //eighth set of harmonics 
  bJSPCDoMixed(kFALSE),
  bJSPCDifferentCharge(kTRUE),
  bJSPCSetSameChargePositive(kTRUE),
  fJSPCMixedHarmonic(0),
  bJSPCComputeEtaGap(kFALSE),		// Do eta gap computation if kTRUE. Default kFALSE
  fJSPCEtaGap(0.8)			// Value of eta gap
{
	// Constructor
	AliInfo("---- AliJSPCTask Constructor ----");
	DefineOutput (1, TList::Class());
}

//____________________________________________________________________________
AliJSPCTask::AliJSPCTask(const AliJSPCTask& ap) :
	AliAnalysisTaskSE(ap.GetName()), 
	fJCatalystTask(ap.fJCatalystTask),
	fJCatalystTaskName(ap.fJCatalystTaskName),
	fIsMC(ap.fIsMC),
	fSPC(ap.fSPC),
  bJSPCSaveAllQA(ap.bJSPCSaveAllQA),
  fJSPCcent_0(ap.fJSPCcent_0), fJSPCcent_1(ap.fJSPCcent_1), fJSPCcent_2(ap.fJSPCcent_2), fJSPCcent_3(ap.fJSPCcent_3), fJSPCcent_4(ap.fJSPCcent_4), fJSPCcent_5(ap.fJSPCcent_5), fJSPCcent_6(ap.fJSPCcent_6), fJSPCcent_7(ap.fJSPCcent_7), fJSPCcent_8(ap.fJSPCcent_8), fJSPCcent_9(ap.fJSPCcent_9), 
  fJSPCcent_10(ap.fJSPCcent_10), fJSPCcent_11(ap.fJSPCcent_11), fJSPCcent_12(ap.fJSPCcent_12), fJSPCcent_13(ap.fJSPCcent_13), fJSPCcent_14(ap.fJSPCcent_14), fJSPCcent_15(ap.fJSPCcent_15), fJSPCcent_16(ap.fJSPCcent_16),
  bJSPCDoFisherYates(ap.bJSPCDoFisherYates),
  fJSPCFisherYatesCutOff(ap.fJSPCFisherYatesCutOff),
  fJSPCMinNumberPart(ap.fJSPCMinNumberPart),
  bJSPCUseWeightsNUE(ap.bJSPCUseWeightsNUE), bJSPCUseWeightsNUA(ap.bJSPCUseWeightsNUA),
  fJSPCNumber(ap.fJSPCNumber),     //number of correlation first correlator
  fJSPCNumberSecond(ap.fJSPCNumberSecond),   //number of correlation second correlator
  fJSPCNumberThird(ap.fJSPCNumberThird),  //number of correlation second correlator
  fJSPCNumberFourth(ap.fJSPCNumberFourth), //number of correlation fourth correlator
  fJSPCNumberFifth(ap.fJSPCNumberFifth),  //number of correlation fifth correlator
  fJSPCNumberSixth(ap.fJSPCNumberSixth),  //number of correlation sixth correlator
  fJSPCNumberSeventh(ap.fJSPCNumberSeventh),  //number of correlation seventh correlator
  fJSPCNumberEighth(ap.fJSPCNumberEighth), //number of correlation eighth correlator
  fJSPCa1(ap.fJSPCa1), fJSPCa2(ap.fJSPCa2), fJSPCa3(ap.fJSPCa3), fJSPCa4(ap.fJSPCa4), fJSPCa5(ap.fJSPCa5), fJSPCa6(ap.fJSPCa6), fJSPCa7(ap.fJSPCa7), //first set of harmonics
  fJSPCb1(ap.fJSPCb1), fJSPCb2(ap.fJSPCb2), fJSPCb3(ap.fJSPCb3), fJSPCb4(ap.fJSPCb4), fJSPCb5(ap.fJSPCb5), fJSPCb6(ap.fJSPCb6), fJSPCb7(ap.fJSPCb7), //second set of harmonics
  fJSPCd1(ap.fJSPCd1), fJSPCd2(ap.fJSPCd2), fJSPCd3(ap.fJSPCd3), fJSPCd4(ap.fJSPCd4), fJSPCd5(ap.fJSPCd5), fJSPCd6(ap.fJSPCd6), fJSPCd7(ap.fJSPCd7), //third set of harmonics
  fJSPCe1(ap.fJSPCe1), fJSPCe2(ap.fJSPCe2), fJSPCe3(ap.fJSPCe3), fJSPCe4(ap.fJSPCe4), fJSPCe5(ap.fJSPCe5), fJSPCe6(ap.fJSPCe6), fJSPCe7(ap.fJSPCe7), //fourth set of harmonics
  fJSPCf1(ap.fJSPCf1), fJSPCf2(ap.fJSPCf2), fJSPCf3(ap.fJSPCf3), fJSPCf4(ap.fJSPCf4), fJSPCf5(ap.fJSPCf5), fJSPCf6(ap.fJSPCf6), fJSPCf7(ap.fJSPCf7),  //fifth set of harmonics 
  fJSPCg1(ap.fJSPCg1), fJSPCg2(ap.fJSPCg2), fJSPCg3(ap.fJSPCg3), fJSPCg4(ap.fJSPCg4), fJSPCg5(ap.fJSPCg5), fJSPCg6(ap.fJSPCg6), fJSPCg7(ap.fJSPCg7),  // sixth set of harmonics 
  fJSPCh1(ap.fJSPCh1), fJSPCh2(ap.fJSPCh2), fJSPCh3(ap.fJSPCh3), fJSPCh4(ap.fJSPCh4), fJSPCh5(ap.fJSPCh5), fJSPCh6(ap.fJSPCh6), fJSPCh7(ap.fJSPCh7), //seventh set of harmonics
  fJSPCi1(ap.fJSPCi1), fJSPCi2(ap.fJSPCi2), fJSPCi3(ap.fJSPCi3), fJSPCi4(ap.fJSPCi4), fJSPCi5(ap.fJSPCi5), fJSPCi6(ap.fJSPCi6), fJSPCi7(ap.fJSPCi7),  //eighth set of harmonics 
  bJSPCDoMixed(ap.bJSPCDoMixed),
  bJSPCDifferentCharge(ap.bJSPCDifferentCharge),
  bJSPCSetSameChargePositive(ap.bJSPCSetSameChargePositive),
  fJSPCMixedHarmonic(ap.fJSPCMixedHarmonic),
  bJSPCComputeEtaGap(ap.bJSPCComputeEtaGap),		// Do eta gap computation if kTRUE. Default kFALSE
  fJSPCEtaGap(ap.fJSPCEtaGap)		// Value of eta gap
{ 

	AliInfo("----DEBUG AliJSPCTask COPY ----");

}

//_____________________________________________________________________________
AliJSPCTask& AliJSPCTask::operator = (const AliJSPCTask& ap)
{
	// assignment operator

	AliInfo("----DEBUG AliJSPCTask operator= ----");
	this->~AliJSPCTask();
	new(this) AliJSPCTask(ap);
	return *this;
}

//______________________________________________________________________________
AliJSPCTask::~AliJSPCTask()
{
	// destructor
	delete fSPC;
}

//________________________________________________________________________

void AliJSPCTask::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fDebug > 1) printf("AliJSPCTask::UserCreateOutPutData() \n");

	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	cout<< "AliJCatalystTask Name = " << fJCatalystTaskName << endl;
	fJCatalystTask = (AliJCatalystTask*)(man->GetTask( fJCatalystTaskName ));

	fSPC = new AliAnalysisSPC("SymmetryPlaneCorrelations",kFALSE);
  fSPC->SetCentrality(fJSPCcent_0, fJSPCcent_1, fJSPCcent_2, fJSPCcent_3, fJSPCcent_4, fJSPCcent_5, fJSPCcent_6, fJSPCcent_7, fJSPCcent_8, fJSPCcent_9, fJSPCcent_10, fJSPCcent_11, fJSPCcent_12, fJSPCcent_13, fJSPCcent_14, fJSPCcent_15, fJSPCcent_16);
  fSPC->SetInitializeCentralityArray();

  fSPC->SetMinNuPar(fJSPCMinNumberPart);
  fSPC->SetFisherYates(bJSPCDoFisherYates, fJSPCFisherYatesCutOff); 

  fSPC->SetUseWeights(bJSPCUseWeightsNUE, bJSPCUseWeightsNUA);

  fSPC->SetCorrSet1(fJSPCNumber, fJSPCa1, fJSPCa2, fJSPCa3, fJSPCa4, fJSPCa5, fJSPCa6, fJSPCa7);
  fSPC->SetCorrSet2(fJSPCNumberSecond, fJSPCb1, fJSPCb2, fJSPCb3, fJSPCb4, fJSPCb5, fJSPCb6, fJSPCb7);
  fSPC->SetCorrSet3(fJSPCNumberThird, fJSPCd1, fJSPCd2, fJSPCd3, fJSPCd4, fJSPCd5, fJSPCd6, fJSPCd7);
  fSPC->SetCorrSet4(fJSPCNumberFourth, fJSPCe1, fJSPCe2, fJSPCe3, fJSPCe4, fJSPCe5, fJSPCe6, fJSPCe7);
  fSPC->SetCorrSet5(fJSPCNumberFifth, fJSPCf1, fJSPCf2, fJSPCf3, fJSPCf4, fJSPCf5, fJSPCf6, fJSPCf7);
  fSPC->SetCorrSet6(fJSPCNumberSixth, fJSPCg1, fJSPCg2, fJSPCg3, fJSPCg4, fJSPCg5, fJSPCg6, fJSPCg7);
  fSPC->SetCorrSet7(fJSPCNumberSeventh, fJSPCh1, fJSPCh2, fJSPCh3, fJSPCh4, fJSPCh5, fJSPCh6, fJSPCh7);
  fSPC->SetCorrSet8(fJSPCNumberEighth, fJSPCi1, fJSPCi2, fJSPCi3, fJSPCi4, fJSPCi5, fJSPCi6, fJSPCi7);
  fSPC->SetMixed(kFALSE,2., kFALSE, kTRUE);
  fSPC->SetEtaGaps(bJSPCComputeEtaGap,fJSPCEtaGap);

	fSPC->SetDebugLevel(fDebug);
	OpenFile(1);

	fSPC->UserCreateOutputObjects();
	PostData(1, fSPC->GetMainList());
}

//______________________________________________________________________________
void AliJSPCTask::UserExec(Option_t* /*option*/) 
{
	// Processing of one event
	if(fDebug > 5) cout << "------- AliJSPCTask Exec-------"<<endl;
	if(!((Entry()-1)%100))  AliInfo(Form(" Processing event # %lld",  Entry()));

	if( fJCatalystTask->GetJCatalystEntry() != fEntry ) return;
	if(fJCatalystTask->GetCentrality()>80. || fJCatalystTask->GetCentrality()<0.) return;
	if(fDebug > 5) cout << Form("Event %d:%d",fEntry, fJCatalystTask->GetJCatalystEntry()) << endl;
	if(fDebug > 5) cout << Form("%s, Nch = %d, cent=%.0f",GetJCatalystTaskName().Data(), fJCatalystTask->GetInputList()->GetEntriesFast(), fJCatalystTask->GetCentrality()) << endl;

	fSPC->SetInputList( fJCatalystTask->GetInputList() );
	fSPC->SetEventCentrality( fJCatalystTask->GetCentrality() );
	fSPC->UserExec("");
}

//______________________________________________________________________________
void AliJSPCTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 

}

//______________________________________________________________________________
void AliJSPCTask::Terminate(Option_t *)
{
	// Processing when the event loop is ended
	cout<<"AliJSPCTask Analysis DONE !!"<<endl; 
}
