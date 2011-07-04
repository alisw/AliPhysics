/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$   */

//  *************************************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  contained in a DB
//  -------------------------------------------------------------
//  W. Ferrarese + P. Cerello Feb 2008
//  INFN Torino
//  SSD QA part: P. Christakoglou - NIKHEF/UU

// --- ROOT system ---
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliITS.h"
#include "AliITSmodule.h"
#include "AliITShit.h"
#include "AliITSdigit.h"
#include "AliITSpListItem.h"
#include "AliRun.h"
#include "AliITSQADataMakerSim.h"
#include "AliITSQASSDDataMakerSim.h"
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"

ClassImp(AliITSQASSDDataMakerSim)

//____________________________________________________________________________ 
AliITSQASSDDataMakerSim::AliITSQASSDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim) :
TObject(),
fAliITSQADataMakerSim(aliITSQADataMakerSim),
//fSSDhTask(0),
fSSDhHTask(0),
fSSDhSTask(0),
fSSDhDTask(0),
fGenOffsetH(0),
fGenOffsetS(0), 
fGenOffsetD(0) 
{
  //ctor used to discriminate OnLine-Offline analysis  
  fGenOffsetH=  new Int_t[AliRecoParam::kNSpecies];                       
  fGenOffsetS=  new Int_t[AliRecoParam::kNSpecies];                           
  fGenOffsetD=  new Int_t[AliRecoParam::kNSpecies];
  for(Int_t i=0; i<AliRecoParam::kNSpecies; i++) 
    {
      fGenOffsetH[i]= 0;
      fGenOffsetS[i]= 0;
      fGenOffsetD[i]= 0;
    }  
}

//____________________________________________________________________________ 
AliITSQASSDDataMakerSim::AliITSQASSDDataMakerSim(const AliITSQASSDDataMakerSim& qadm) :
TObject(),
fAliITSQADataMakerSim(qadm.fAliITSQADataMakerSim),
//fSSDhTask(qadm.fSSDhTask),
fSSDhHTask(qadm.fSSDhHTask),
fSSDhSTask(qadm.fSSDhSTask),
fSSDhDTask(qadm.fSSDhDTask),
fGenOffsetH(qadm.fGenOffsetH), 
fGenOffsetS(qadm.fGenOffsetS), 
fGenOffsetD(qadm.fGenOffsetD) 
{
  //copy ctor 
  fAliITSQADataMakerSim->SetName((const char*)qadm.fAliITSQADataMakerSim->GetName()) ; 
  fAliITSQADataMakerSim->SetTitle((const char*)qadm.fAliITSQADataMakerSim->GetTitle());
  }

//__________________________________________________________________
AliITSQASSDDataMakerSim& AliITSQASSDDataMakerSim::operator = (const AliITSQASSDDataMakerSim& qac ) {
  // Equal operator.
  this->~AliITSQASSDDataMakerSim();
  new(this) AliITSQASSDDataMakerSim(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerSim::StartOfDetectorCycle() {
  //Detector specific actions at start of cycle
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Start of SSD Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t /*task*/, TObjArray** /*list*/) {
  // launch the QA checking
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSDM instantiates checker with Run(AliQAv1::kITS, task, list)\n"); 
  
//  AliQAChecker::Instance()->Run( AliQAv1::kITS , task, list);
}

//____________________________________________________________________________ 
Int_t AliITSQASSDDataMakerSim::InitDigits() { 
  // Initialization for DIGIT data - SSD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Int_t rv = 0 ; 
 // fGenOffsetD = (fAliITSQADataMakerSim->fDigitsQAList[AliRecoParam::kDefault])->GetEntries();

  // custom code here
  TH1F *fHistSSDModule = new TH1F("fHistSSDDigitsModule",
				  "SSD Digits Module;SSD Module Number;N_{DIGITS}",
				  1698,499.5,2197.5);  
  rv = fAliITSQADataMakerSim->Add2DigitsList(fHistSSDModule,
					fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()] + 0, !expert, image);
  fSSDhDTask += 1;
  TH2F *fHistSSDModuleStrip = new TH2F("fHistSSDDigitsModuleStrip",
				       "SSD Digits Module Strip;N_{Strip};N_{Module}",
				       1540,0,1540,1698,499.5,2197.5);  
  rv = fAliITSQADataMakerSim->Add2DigitsList(fHistSSDModuleStrip,
					fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()] + 1, !expert, image);
  fSSDhDTask += 1;

  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SSD Digits histograms booked\n",fSSDhDTask));
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASSDDataMakerSim::MakeDigits(TTree *digits) 
{ 
  // Fill QA for DIGIT - SSD -
  Int_t rv = 0 ; 

  AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  fITS->SetTreeAddress();
  TClonesArray *iSSDdigits  = fITS->DigitsAddress(2);
  for(Int_t iModule = 500; iModule < 2198; iModule++) {
    iSSDdigits->Clear();
    digits->GetEvent(iModule);    
    Int_t ndigits = iSSDdigits->GetEntries();
    fAliITSQADataMakerSim->FillDigitsData(fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()] + 0,iModule,ndigits);
    if(ndigits != 0)
      AliDebug(AliQAv1::GetQADebugLevel(),Form("Module: %d - Digits: %d",iModule,ndigits));
 
    for (Int_t iDigit = 0; iDigit < ndigits; iDigit++) {
      AliITSdigit *dig = (AliITSdigit*)iSSDdigits->UncheckedAt(iDigit);
      Int_t fStripNumber = (dig->GetCoord1() == 0) ? dig->GetCoord2() : dig->GetCoord2() + fgkNumberOfPSideStrips;
      fAliITSQADataMakerSim->FillDigitsData(fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()] + 1,fStripNumber,iModule,dig->GetSignal());
    }//digit loop
  }//module loop
  //
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASSDDataMakerSim::InitSDigits() { 
  // Initialization for SDIGIT data - SSD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Int_t rv = 0 ; 
  //fGenOffsetS = (fAliITSQADataMakerSim->fSDigitsQAList[AliRecoParam::kDefault])->GetEntries();

  // custom code here
  TH1F *fHistSSDModule = new TH1F("fHistSSDSDigitsModule",
				  "SSD SDigits Module;SSD Module Number;N_{SDIGITS}",
				  1698,499.5,2197.5);  
  rv = fAliITSQADataMakerSim->Add2SDigitsList(fHistSSDModule,
					 fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()] + 0, !expert, image);
  fSSDhSTask += 1;  

  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SSD SDigits histograms booked\n",fSSDhSTask));
  //
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASSDDataMakerSim::MakeSDigits(TTree *sdigits) { 
  // Fill QA for SDIGIT - SSD -
  Int_t rv = 0 ; 
 
  static TClonesArray iSSDEmpty("AliITSpListItem",10000);
  iSSDEmpty.Clear();
  TClonesArray *iSSDsdigits = &iSSDEmpty;

  //  AliDebug(AliQAv1::GetQADebugLevel(), Form("Trying to access the sdigits histogram: %d\n",fGenOffsetS));

  TBranch *brchSDigits = sdigits->GetBranch("ITS");
  brchSDigits->SetAddress(&iSSDsdigits);
  for(Int_t iModule = 500; iModule < 2198; iModule++) {
    iSSDsdigits->Clear();
    sdigits->GetEvent(iModule);    
    Int_t ndigits = iSSDsdigits->GetEntries();
    fAliITSQADataMakerSim->FillSDigitsData(fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()] + 0,iModule,ndigits);
    if(ndigits != 0)
      AliDebug(AliQAv1::GetQADebugLevel(),Form("Module: %d - Digits: %d",iModule,ndigits));

//     for (Int_t iDigit = 0; iDigit < ndigits; iDigit++) {
//       AliITSpListItem *dig=(AliITSpListItem*)iSSDsdigits->At(iDigit);
//       dig=0;
//     }//digit loop
  }//module loop
  //
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASSDDataMakerSim::InitHits() { 
  // Initialization for HITS data - SSD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Int_t rv = 0 ; 

  //fGenOffsetH = (fAliITSQADataMakerSim->fHitsQAList[fEventSpecie])->GetEntries();

  // custom code here
  TH1F *fHistSSDModule = new TH1F("fHistSSDHitsModule",
				  "SSD Hits Module;SSD Module Number;N_{HITS}",
				  1698,499.5,2197.5); 
  rv = fAliITSQADataMakerSim->Add2HitsList(fHistSSDModule,
				      fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 0, !expert, image);
  fSSDhHTask += 1;
  TH1F *fHistSSDGlobalX = new TH1F("fHistSSDHitsGlobalX",
				   "SSD Hits Global X;x [cm];Entries",
				   1000,-50.,50.);
  rv = fAliITSQADataMakerSim->Add2HitsList(fHistSSDGlobalX,
				      fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 1, !expert, image);
  fSSDhHTask += 1;
  TH1F *fHistSSDGlobalY = new TH1F("fHistSSDHitsGlobalY",
				   "SSD Hits Global Y;y [cm];Entries",
				   1000,-50.,50.);
  rv = fAliITSQADataMakerSim->Add2HitsList(fHistSSDGlobalY,
				      fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 2, !expert, image);
  fSSDhHTask += 1;
  TH1F *fHistSSDGlobalZ = new TH1F("fHistSSDHitsGlobalZ",
				   "SSD Hits Global Z ;z [cm];Entries",
				   1000,-60.,60.);
  rv = fAliITSQADataMakerSim->Add2HitsList(fHistSSDGlobalZ,
				      fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 3, !expert, image);
  fSSDhHTask += 1;
  TH1F *fHistSSDLocalX = new TH1F("fHistSSDHitsLocalX",
				  "SSD Hits Local X;x [cm];Entries",
				  1000,-4.,4.);
  rv = fAliITSQADataMakerSim->Add2HitsList(fHistSSDLocalX,
				      fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 4, !expert, image);
  fSSDhHTask += 1;
  TH1F *fHistSSDLocalY = new TH1F("fHistSSDHitsLocalY",
				  "SSD Hits Local Y;y [cm];Entries",
				  1000,-0.1,0.1);
  rv = fAliITSQADataMakerSim->Add2HitsList(fHistSSDLocalY,
				      fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 5, !expert, image);
  fSSDhHTask += 1;
  TH1F *fHistSSDLocalZ = new TH1F("fHistSSDHitsLocalZ",
				  "SSD Hits Local Z;z [cm];Entries",
				  1000,-4.,4.);
  rv = fAliITSQADataMakerSim->Add2HitsList(fHistSSDLocalZ,
				      fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 6, !expert, image);
  fSSDhHTask += 1;
  TH1F *fHistSSDIonization = new TH1F("fHistSSDHitsIonization",
				      "SSD Hits Ionization;log(dE/dx) [KeV];N_{Hits}",
				      100,-7,-2);
  rv = fAliITSQADataMakerSim->Add2HitsList(fHistSSDIonization,
				      fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 7, !expert, image);
  fSSDhHTask += 1;
  TH2F *fHistSSDGlobalXY = new TH2F("fHistSSDHitsGlobalXY",
				    "SSD Hits Global XY;x [cm];y [cm]",
				    1000,-50.,50.,
				    1000,-50.,50.);
  rv = fAliITSQADataMakerSim->Add2HitsList(fHistSSDGlobalXY,
				      fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 8, !expert, image);
  fSSDhHTask += 1;
 
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SSD Hits histograms booked\n",fSSDhHTask));
  return rv ; 
}


//____________________________________________________________________________
Int_t AliITSQASSDDataMakerSim::MakeHits(TTree *hits) { 
  // Fill QA for HITS - SSD -
  Int_t rv = 0 ; 
 
  AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  fITS->SetTreeAddress();
  Int_t nmodules;
  fITS->InitModules(-1,nmodules);
  fITS->FillModules(hits,0);
  for(Int_t iModule = 500; iModule < 2198; iModule++) {
    AliITSmodule *module = fITS->GetModule(iModule);
    TObjArray *arrHits = module->GetHits();
    Int_t nhits = arrHits->GetEntriesFast();
    if(nhits != 0)
      AliDebug(AliQAv1::GetQADebugLevel(),Form("Module: %d - Hits: %d",iModule,nhits));
    for (Int_t iHit = 0; iHit < nhits; iHit++) {
      AliITShit *hit = (AliITShit*) arrHits->At(iHit);
      
      fAliITSQADataMakerSim->FillHitsData(fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 0,iModule);
      fAliITSQADataMakerSim->FillHitsData(fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 1,hit->GetXG());
      fAliITSQADataMakerSim->FillHitsData(fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 2,hit->GetYG());
      fAliITSQADataMakerSim->FillHitsData(fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 3,hit->GetZG());
      fAliITSQADataMakerSim->FillHitsData(fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 4,hit->GetXL());
      fAliITSQADataMakerSim->FillHitsData(fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 5,hit->GetYL());
      fAliITSQADataMakerSim->FillHitsData(fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 6,hit->GetZL());
      if(hit->GetIonization())
	fAliITSQADataMakerSim->FillHitsData(fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 7,TMath::Log10(hit->GetIonization()));
      fAliITSQADataMakerSim->FillHitsData(fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()] + 8,hit->GetXG(),hit->GetYG());
    }//hit loop
  }//module loop  
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASSDDataMakerSim::GetOffset(AliQAv1::TASKINDEX_t task,Int_t specie){
  // Returns histogram offset according to the specified task
  Int_t offset=0;
  if( task == AliQAv1::kHITS){
    offset=fGenOffsetH[specie];  
  }
  else if( task == AliQAv1::kSDIGITS) {
    offset=fGenOffsetS[specie];   
  }
  else if( task == AliQAv1::kDIGITS) {
    offset=fGenOffsetD[specie];   
  }
  else {
    AliInfo("No task has been selected. TaskHisto set to zero.\n");
  }

  return offset;
}


//____________________________________________________________________________ 
void AliITSQASSDDataMakerSim::SetOffset(AliQAv1::TASKINDEX_t task, Int_t offset,Int_t specie ){
  // Returns histogram offset according to the specified task
  if( task == AliQAv1::kHITS){
    fGenOffsetH[specie] = offset;  
  }
  else if( task == AliQAv1::kSDIGITS) {
    fGenOffsetS[specie] = offset;   
  }
  else if( task == AliQAv1::kDIGITS) {
    fGenOffsetD[specie] = offset;   
  }
  else {
    AliInfo("No task has been selected. TaskHisto set to zero.\n");
  }
}

//____________________________________________________________________________ 
Int_t AliITSQASSDDataMakerSim::GetTaskHisto(AliQAv1::TASKINDEX_t task) {
  // Returns the number of booked histograms for the selected task
  Int_t histotot=0;
  if( task == AliQAv1::kHITS) {
    histotot=fSSDhHTask ;  
  }
  else if( task == AliQAv1::kSDIGITS) {
    histotot=fSSDhSTask;   
  }
  else if( task == AliQAv1::kDIGITS) {
    histotot=fSSDhDTask ;   
  }
  else {
    AliInfo("No task has been selected. TaskHisto set to zero.\n");
  }
  return histotot;

}
