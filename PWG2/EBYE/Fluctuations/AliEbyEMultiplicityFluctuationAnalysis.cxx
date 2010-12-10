/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Satyajit Jena      .                                           *
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

/*-------------------------------------------------------------------------
 *
 *           AliEbyEMultiplicityFluctuationAnalysis base class
 *        This class deals with the multiplicity  fluctuation 
 *              Origin: Satyajit Jena <sjena@cern.ch>
 *
 ------------------------------------------------------------------------*/

#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH3F.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TParticle.h>
#include <TList.h>

#include <AliExternalTrackParam.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliMCParticle.h>
#include <AliPID.h>
#include <AliStack.h>
#include <AliCFContainer.h>
#include <AliCFEffGrid.h>
#include <AliCFDataGrid.h>
#include <AliTPCPIDResponse.h>
#include <AliESDpid.h>
#include "AliEbyEEventBase.h"
#include "AliEbyEMultiplicityFluctuationAnalysis.h"

ClassImp(AliEbyEMultiplicityFluctuationAnalysis)

//____________________________________________________________________//
AliEbyEMultiplicityFluctuationAnalysis::AliEbyEMultiplicityFluctuationAnalysis() : 
  TObject(), 
  fListMFQA(0),
  fListMeasureMF(0),
  fEbyEBase(0)
{
  //Default constructor
  InitHistos();

}

//____________________________________________________________________//
AliEbyEMultiplicityFluctuationAnalysis::~AliEbyEMultiplicityFluctuationAnalysis() {

  //Default destructor
  if(fEbyEBase) delete fEbyEBase;
  if(fListMFQA) delete fListMFQA;
  if(fListMeasureMF) delete fListMeasureMF;

}

//____________________________________________________________________//
void AliEbyEMultiplicityFluctuationAnalysis::InitHistos() {

  fListMFQA = new TList();
  fListMFQA->SetName("MFQaList");

  fListMeasureMF = new TList();
  fListMeasureMF->SetName("MFMeasureList");

  TH1F *fEvtStat = new TH1F("hEventCounterMH"," Event Counters for MF Analysis", 20,0,20);
  fListMFQA->Add(fEvtStat); // --:0
 
  TH2F *fhMultiplicity = new TH2F("hMultiplicity"," Multiplisity in Centrality Bin", 50,0,50, 3000,0,3000);
  fListMeasureMF->Add(fhMultiplicity); //! -->:0


}


//____________________________________________________________________//
void AliEbyEMultiplicityFluctuationAnalysis::Analyze(AliESDEvent* esd) {

  // Analysis from ESD
  Int_t cent = 0;

  ((TH1F*)(fListMFQA->At(0)))->Fill(1.);  
  Int_t itr = esd->GetNumberOfTracks();
  ((TH2F*)(fListMeasureMF->At(0)))->Fill(cent,itr); 

}


//____________________________________________________________________//
void AliEbyEMultiplicityFluctuationAnalysis::Analyze(AliAODEvent* aod) {
  // Analysis from AOD
  ((TH1F*)(fListMFQA->At(0)))->Fill(1.);  
  printf("%d\n",aod->GetNTracks());


}

//____________________________________________________________________//
void AliEbyEMultiplicityFluctuationAnalysis::Analyze(AliStack* stack) {
  // Analysis from MC stack
  ((TH1F*)(fListMFQA->At(0)))->Fill(1.);

  printf("%d \n",stack->GetNtrack());

  
}
