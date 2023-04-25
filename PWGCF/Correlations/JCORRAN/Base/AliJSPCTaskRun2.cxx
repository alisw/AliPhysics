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
#include "AliJSPCTaskRun2.h"

//______________________________________________________________________________
AliJSPCTaskRun2::AliJSPCTaskRun2() :   
  AliAnalysisTaskSE("JSPCTask"),
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fIsMC(kFALSE),
	fSPC(NULL),
	bJSPCRun2SaveAllQA(kTRUE),
	fAliSPCRun2cent_0(0.), fAliSPCRun2cent_1(0.), fAliSPCRun2cent_2(0.), fAliSPCRun2cent_3(0.), fAliSPCRun2cent_4(0.), fAliSPCRun2cent_5(0.), fAliSPCRun2cent_6(0.), fAliSPCRun2cent_7(0.), fAliSPCRun2cent_8(0.), fAliSPCRun2cent_9(0.), 
  fAliSPCRun2cent_10(0.), fAliSPCRun2cent_11(0.), fAliSPCRun2cent_12(0.), fAliSPCRun2cent_13(0.), fAliSPCRun2cent_14(0.), fAliSPCRun2cent_15(0.), fAliSPCRun2cent_16(0.),
  fAliSPCRun2MinNumberPart(14),
  bAliSPCRun2UseWeightsNUE(kTRUE),
  bAliSPCRun2UseWeightsNUA(kFALSE), 
  bAliSPCRun2ComputeEtaGap(kFALSE),		// Do eta gap computation if kTRUE. Default kFALSE
  fAliSPCRun2EtaGap(0.8)			// Value of eta gap
{
}

//______________________________________________________________________________
AliJSPCTaskRun2::AliJSPCTaskRun2(const char *name):
	AliAnalysisTaskSE(name), 
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fIsMC(kFALSE),
	fSPC(NULL),
	bJSPCRun2SaveAllQA(kTRUE),
	fAliSPCRun2cent_0(0.), fAliSPCRun2cent_1(0.), fAliSPCRun2cent_2(0.), fAliSPCRun2cent_3(0.), fAliSPCRun2cent_4(0.), fAliSPCRun2cent_5(0.), fAliSPCRun2cent_6(0.), fAliSPCRun2cent_7(0.), fAliSPCRun2cent_8(0.), fAliSPCRun2cent_9(0.), 
  fAliSPCRun2cent_10(0.), fAliSPCRun2cent_11(0.), fAliSPCRun2cent_12(0.), fAliSPCRun2cent_13(0.), fAliSPCRun2cent_14(0.), fAliSPCRun2cent_15(0.), fAliSPCRun2cent_16(0.),
  fAliSPCRun2MinNumberPart(14),
  bAliSPCRun2UseWeightsNUE(kTRUE),
  bAliSPCRun2UseWeightsNUA(kFALSE),
  bAliSPCRun2ComputeEtaGap(kFALSE),		// Do eta gap computation if kTRUE. Default kFALSE
  fAliSPCRun2EtaGap(0.8)			// Value of eta gap
{
	// Constructor
	AliInfo("---- AliJSPCTaskRun2 Constructor ----");
	DefineOutput (1, TList::Class());
}

//____________________________________________________________________________
AliJSPCTaskRun2::AliJSPCTaskRun2(const AliJSPCTaskRun2& ap) :
	AliAnalysisTaskSE(ap.GetName()), 
	fJCatalystTask(ap.fJCatalystTask),
	fJCatalystTaskName(ap.fJCatalystTaskName),
	fIsMC(ap.fIsMC),
	fSPC(ap.fSPC),
  bJSPCRun2SaveAllQA(ap.bJSPCRun2SaveAllQA),
  fAliSPCRun2cent_0(ap.fAliSPCRun2cent_0), fAliSPCRun2cent_1(ap.fAliSPCRun2cent_1), fAliSPCRun2cent_2(ap.fAliSPCRun2cent_2), fAliSPCRun2cent_3(ap.fAliSPCRun2cent_3), fAliSPCRun2cent_4(ap.fAliSPCRun2cent_4), fAliSPCRun2cent_5(ap.fAliSPCRun2cent_5), fAliSPCRun2cent_6(ap.fAliSPCRun2cent_6), fAliSPCRun2cent_7(ap.fAliSPCRun2cent_7), fAliSPCRun2cent_8(ap.fAliSPCRun2cent_8), fAliSPCRun2cent_9(ap.fAliSPCRun2cent_9), 
  fAliSPCRun2cent_10(ap.fAliSPCRun2cent_10), fAliSPCRun2cent_11(ap.fAliSPCRun2cent_11), fAliSPCRun2cent_12(ap.fAliSPCRun2cent_12), fAliSPCRun2cent_13(ap.fAliSPCRun2cent_13), fAliSPCRun2cent_14(ap.fAliSPCRun2cent_14), fAliSPCRun2cent_15(ap.fAliSPCRun2cent_15), fAliSPCRun2cent_16(ap.fAliSPCRun2cent_16),
  fAliSPCRun2MinNumberPart(ap.fAliSPCRun2MinNumberPart),
  bAliSPCRun2UseWeightsNUE(ap.bAliSPCRun2UseWeightsNUE), bAliSPCRun2UseWeightsNUA(ap.bAliSPCRun2UseWeightsNUA),
  bAliSPCRun2ComputeEtaGap(ap.bAliSPCRun2ComputeEtaGap),		// Do eta gap computation if kTRUE. Default kFALSE
  fAliSPCRun2EtaGap(ap.fAliSPCRun2EtaGap)		// Value of eta gap
{ 

	AliInfo("----DEBUG AliJSPCTaskRun2 COPY ----");

}

//_____________________________________________________________________________
AliJSPCTaskRun2& AliJSPCTaskRun2::operator = (const AliJSPCTaskRun2& ap)
{
	// assignment operator

	AliInfo("----DEBUG AliJSPCTaskRun2 operator= ----");
	this->~AliJSPCTaskRun2();
	new(this) AliJSPCTaskRun2(ap);
	return *this;
}

//______________________________________________________________________________
AliJSPCTaskRun2::~AliJSPCTaskRun2()
{
	// destructor
	delete fSPC;
}

//________________________________________________________________________

void AliJSPCTaskRun2::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fDebug > 1) printf("AliJSPCTaskRun2::UserCreateOutPutData() \n");

	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	cout<< "AliJCatalystTask Name = " << fJCatalystTaskName << endl;
	fJCatalystTask = (AliJCatalystTask*)(man->GetTask( fJCatalystTaskName ));

	fSPC = new AliAnalysisSPCRun2("SymmetryPlaneCorrelations");
  fSPC->SetCentrality(fAliSPCRun2cent_0, fAliSPCRun2cent_1, fAliSPCRun2cent_2, fAliSPCRun2cent_3, fAliSPCRun2cent_4, fAliSPCRun2cent_5, fAliSPCRun2cent_6, fAliSPCRun2cent_7, fAliSPCRun2cent_8, fAliSPCRun2cent_9, fAliSPCRun2cent_10, fAliSPCRun2cent_11, fAliSPCRun2cent_12, fAliSPCRun2cent_13, fAliSPCRun2cent_14, fAliSPCRun2cent_15, fAliSPCRun2cent_16);
  fSPC->SetInitializeCentralityArray();

  fSPC->SetMinNuPar(fAliSPCRun2MinNumberPart);

  fSPC->SetUseWeights(bAliSPCRun2UseWeightsNUE, bAliSPCRun2UseWeightsNUA);	// LOKI: Decommented

  for (int i=0; i<12; i++) {
    fSPC->SetCorrSet(i, fAliSPCRun2HarmosArray[i]);
  }

  fSPC->SetEtaGaps(bAliSPCRun2ComputeEtaGap,fAliSPCRun2EtaGap);

	fSPC->SetDebugLevel(fDebug);
	OpenFile(1);

	fSPC->UserCreateOutputObjects();
	PostData(1, fSPC->GetMainList());
}

//______________________________________________________________________________
void AliJSPCTaskRun2::UserExec(Option_t* /*option*/) 
{
	// Processing of one event
	if(fDebug > 5) cout << "------- AliJSPCTaskRun2 Exec-------"<<endl;
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
void AliJSPCTaskRun2::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 

}

//______________________________________________________________________________
void AliJSPCTaskRun2::Terminate(Option_t *)
{
	// Processing when the event loop is ended
	cout<<"AliJSPCTaskRun2 Analysis DONE !!"<<endl; 
}
