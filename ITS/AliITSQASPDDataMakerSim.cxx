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

/* $Id$  */
//  *************************************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  contained in a DB
//  -------------------------------------------------------------
//  W. Ferrarese + P. Cerello INFN Torino Feb 2008
//  M. Nicassio D. Elia INFN Bari April 2008
//  maria.nicassio@ba.infn.it

// --- ROOT system ---
#include <TTree.h>
#include <TH2.h>
#include <TH1.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliITSQADataMakerSim.h"
#include "AliITSQASPDDataMakerSim.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliITSdigit.h"   
#include "AliITSdigitSPD.h"
#include "AliITS.h"
#include "AliITSmodule.h"
#include "AliITShit.h"
#include "AliITSLoader.h"
#include "AliRunLoader.h"

ClassImp(AliITSQASPDDataMakerSim)

//____________________________________________________________________________ 
AliITSQASPDDataMakerSim::AliITSQASPDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim) :
TObject(),
fAliITSQADataMakerSim(aliITSQADataMakerSim),
fSPDhHTask(0),
fSPDhSTask(0),
fSPDhDTask(0),
fGenOffsetH(0),
fGenOffsetS(0),
fGenOffsetD(0)
{
  //ctor used to discriminate OnLine-Offline analysis   
  fGenOffsetH=  new Int_t[AliRecoParam::kNSpecies];                       
  fGenOffsetS=  new Int_t[AliRecoParam::kNSpecies];                           
  fGenOffsetD=  new Int_t[AliRecoParam::kNSpecies];
  for(Int_t i=0; i<AliRecoParam::kNSpecies; i++) {
    fGenOffsetH[i]= 0;
    fGenOffsetS[i]= 0;
    fGenOffsetD[i]= 0;
  }             
}

//____________________________________________________________________________ 
AliITSQASPDDataMakerSim::AliITSQASPDDataMakerSim(const AliITSQASPDDataMakerSim& qadm) :
TObject(),
fAliITSQADataMakerSim(qadm.fAliITSQADataMakerSim),
fSPDhHTask(qadm.fSPDhHTask),
fSPDhSTask(qadm.fSPDhSTask),
fSPDhDTask(qadm.fSPDhDTask),
fGenOffsetH(qadm.fGenOffsetH),
fGenOffsetS(qadm.fGenOffsetS),
fGenOffsetD(qadm.fGenOffsetD)
{
  //copy ctor 
  fAliITSQADataMakerSim->SetName((const char*)qadm.fAliITSQADataMakerSim->GetName()) ; 
  fAliITSQADataMakerSim->SetTitle((const char*)qadm.fAliITSQADataMakerSim->GetTitle());

}

//__________________________________________________________________
AliITSQASPDDataMakerSim& AliITSQASPDDataMakerSim::operator = (const AliITSQASPDDataMakerSim& qac )
{
  // Equal operator.
  this->~AliITSQASPDDataMakerSim();
  new(this) AliITSQASPDDataMakerSim(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Start of SPD Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSDM instantiates checker with Run(AliQAv1::kITS, task, list)\n"); 
  
  //AliQAChecker::Instance()->Run( AliQAv1::kITS , task, list);
}

//____________________________________________________________________________ 
Int_t AliITSQASPDDataMakerSim::InitDigits()
{ 
  // Initialization for DIGIT data - SPD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ;
  Int_t rv = 0 ; 
  //fGenOffsetD = (fAliITSQADataMakerSim->fDigitsQAList[AliRecoParam::kDefault])->GetEntries();
  //fSPDhDTask must be incremented by one unit every time a histogram is ADDED to the QA List

  Char_t name[50];
  Char_t title[50];

  TH1F *hlayer = new TH1F("SPDLayPattern_SPD","Layer map - SPD",6,0.,6.);
  hlayer->GetXaxis()->SetTitle("Layer number");
  hlayer->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerSim->Add2DigitsList(hlayer,fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()], expert, !image);
  fSPDhDTask++;
  
  TH1F **hmod = new TH1F*[2];
  for (Int_t iLay=0; iLay<2; iLay++) {
    sprintf(name,"SPDModPattern_SPD%d",iLay+1);
    sprintf(title,"Module map - SPD Layer %d",iLay+1);
    hmod[iLay]=new TH1F(name,title,240,0,240);
    hmod[iLay]->GetXaxis()->SetTitle("Module number");
    hmod[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerSim->Add2DigitsList(hmod[iLay],1+iLay+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
    fSPDhDTask++;
  }
  
  TH1F *hcolumns = new TH1F("SPDColumns_SPD","Columns - SPD",160,0.,160.);
  hcolumns->GetXaxis()->SetTitle("Column number");
  hcolumns->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerSim->Add2DigitsList(hcolumns,3+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()], expert, !image);
  fSPDhDTask++;

  TH1F *hrows = new TH1F("SPDRows_SPD","Rows - SPD",256,0.,256.);
  hrows->GetXaxis()->SetTitle("Row number");
  hrows->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerSim->Add2DigitsList(hrows,4+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()], expert, !image);
  fSPDhDTask++;

  TH1F** hMultSPDdigits = new TH1F*[2];
  for (Int_t iLay=0; iLay<2; ++iLay) {
    sprintf(name,"SPDDigitMultiplicity_SPD%d",iLay+1);
    sprintf(title,"Digit multiplicity - SPD Layer %d",iLay+1);
    hMultSPDdigits[iLay]=new TH1F(name,title,200,0.,200.);
    hMultSPDdigits[iLay]->GetXaxis()->SetTitle("Digit multiplicity");
    hMultSPDdigits[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerSim->Add2DigitsList(hMultSPDdigits[iLay], 5+iLay+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
    fSPDhDTask++;
  }

  TH2F *hMultSPDdig2MultSPDdig1 
       = new TH2F("SPDDigitMultCorrelation_SPD","Digit multiplicity correlation - SPD",200,0.,200.,200,0.,200.);
  hMultSPDdig2MultSPDdig1->GetXaxis()->SetTitle("Digit multiplicity (Layer 1)");
  hMultSPDdig2MultSPDdig1->GetYaxis()->SetTitle("Digit multiplicity (Layer 2)");
  rv = fAliITSQADataMakerSim->Add2DigitsList(hMultSPDdig2MultSPDdig1,7+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSPDhDTask++;

  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SPD Digits histograms booked\n",fSPDhDTask));
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASPDDataMakerSim::MakeDigits(TTree *digits)
{ 
  // Fill QA for DIGIT - SPD -
  Int_t rv = 0 ; 
 
  AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  fITS->SetTreeAddress();
  TClonesArray *iITSdigits  = fITS->DigitsAddress(0);  // 0->SPD

  Int_t nDigitsL1=0;
  Int_t nDigitsL2=0;

  for (Int_t imod=0; imod<240; ++imod){
    digits->GetEvent(imod);
    Int_t ndigits = iITSdigits->GetEntries();
    if (imod<80) {
      fAliITSQADataMakerSim->GetDigitsData(0+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(0.5,ndigits);
      fAliITSQADataMakerSim->GetDigitsData(1+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(imod,ndigits);
      nDigitsL1+=ndigits;
    }
    else {
      fAliITSQADataMakerSim->GetDigitsData(0+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(1,ndigits);
      fAliITSQADataMakerSim->GetDigitsData(2+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(imod,ndigits);
      nDigitsL2+=ndigits;
    }
    for (Int_t idig=0; idig<ndigits; ++idig) {
      AliITSdigit *dig=(AliITSdigit*)iITSdigits->UncheckedAt(idig);
      Int_t col=dig->GetCoord1();  // cell number z
      Int_t row=dig->GetCoord2();  // cell number x
      fAliITSQADataMakerSim->GetDigitsData(3+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(col);
      fAliITSQADataMakerSim->GetDigitsData(4+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(row);
    }
  }
  fAliITSQADataMakerSim->GetDigitsData(5+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(nDigitsL1);
  fAliITSQADataMakerSim->GetDigitsData(6+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(nDigitsL2);
  fAliITSQADataMakerSim->GetDigitsData(7+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(nDigitsL1,nDigitsL2);
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASPDDataMakerSim::InitSDigits()
{ 
  // Initialization for SDIGIT data - SPD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ;
  Int_t rv = 0 ; 
  //fGenOffsetS = (fAliITSQADataMakerSim->fSDigitsQAList[AliRecoParam::kDefault])->GetEntries();
  //printf("--W-- AliITSQASPDDataMakerSim::InitSDigits()  fGenOffset= %d \n",fGenOffset);
  //fSPDhSTask must be incremented by one unit every time a histogram is ADDED to the QA List
  
  Char_t name[50];
  Char_t title[50];

  TH1F *hlayer = new TH1F("SPDLayPattern_SPD","Layer map - SPD",6,0.,6.);
  hlayer->GetXaxis()->SetTitle("Layer number");
  hlayer->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerSim->Add2SDigitsList(hlayer,fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()], expert, !image);
  fSPDhSTask++;

  TH1F **hmod = new TH1F*[2];
  for (Int_t iLay=0; iLay<2; ++iLay) {
    sprintf(name,"SPDModPattern_SPD%d",iLay+1);
    sprintf(title,"Module map - SPD Layer %d",iLay+1);
    hmod[iLay]=new TH1F(name,title,240,0,240);
    hmod[iLay]->GetXaxis()->SetTitle("Module number");
    hmod[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerSim->Add2SDigitsList(hmod[iLay],1+iLay+fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
    fSPDhSTask++;
  }

  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SPD SDigits histograms booked\n",fSPDhSTask));
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASPDDataMakerSim::MakeSDigits(TTree *sdigits)
{ 
  // Fill QA for SDIGIT - SPD -
  Int_t rv = 0 ; 
  
  static TClonesArray * sdig ; 
  if (! sdig )
    sdig = new TClonesArray( "AliITSpListItem",1000 );
  
  TBranch *brchSDigits = sdigits->GetBranch("ITS");
  for (Int_t imod=0; imod<240; ++imod){
    brchSDigits->SetAddress( &sdig );
    brchSDigits->GetEvent(imod);
    Int_t nsdig=sdig->GetEntries();
    if (imod<80) {
      fAliITSQADataMakerSim->GetSDigitsData(0+fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(0.5,nsdig);
      fAliITSQADataMakerSim->GetSDigitsData(1+fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(imod,nsdig);
    }
    else {
      fAliITSQADataMakerSim->GetSDigitsData(0+fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(1,nsdig);
      fAliITSQADataMakerSim->GetSDigitsData(2+fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(imod,nsdig);
    }
    sdig->Clear() ;
  }
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASPDDataMakerSim::InitHits()
{ 
  // Initialization for HITS data - SPD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ;
  Int_t rv = 0 ; 
  
  //fGenOffsetH = (fAliITSQADataMakerSim->fHitsQAList[AliRecoParam::kDefault])->GetEntries();
  //printf("--W-- AliITSQASPDDataMakerSim::InitHits()  fGenOffset= %d \n",fGenOffset);
  //fSPDhHTask must be incremented by one unit every time a histogram is ADDED to the QA List
  Char_t name[50];
  Char_t title[50];

  TH1F *hlayer = new TH1F("SPDLayPattern_SPD","Layer map - SPD",6,0.,6.);
  hlayer->GetXaxis()->SetTitle("Layer number");
  hlayer->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerSim->Add2HitsList(hlayer,fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()], expert, !image);
  fSPDhHTask++;

  TH1F **hmod = new TH1F*[2];
  for (Int_t iLay=0; iLay<2; ++iLay) {
    sprintf(name,"SPDModPattern_SPD%d",iLay+1);
    sprintf(title,"Module map - SPD Layer %d",iLay+1);
    hmod[iLay]=new TH1F(name,title,240,0,240);
    hmod[iLay]->GetXaxis()->SetTitle("Module number");
    hmod[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerSim->Add2HitsList(hmod[iLay],1+iLay+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
    fSPDhHTask++;
  }

  TH1F *hhitlenght = new TH1F("SPDLenght_SPD","SPD Hit lenght along y_{loc} coord",210,0.,210.);
  hhitlenght->GetXaxis()->SetTitle("Hit lenght [#mum]");
  hhitlenght->GetYaxis()->SetTitle("# hits");
  rv = fAliITSQADataMakerSim->Add2HitsList(hhitlenght,3+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSPDhHTask++;

  TH1F *hEdepos = new TH1F("SPDEnergyDeposit_SPD","SPD Deposited energy distribution (y_{loc}>180 #mum)",150,0.,300.); 
  hEdepos->GetXaxis()->SetTitle("Deposited energy [keV]"); 
  hEdepos->GetYaxis()->SetTitle("# hits");
  rv = fAliITSQADataMakerSim->Add2HitsList(hEdepos,4+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSPDhHTask++;

  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SPD Hits histograms booked\n",fSPDhHTask));
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASPDDataMakerSim::MakeHits(TTree *hits)
{ 
  // Fill QA for HITS - SPD -
  Int_t rv = 0 ; 
 
  AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  fITS->SetTreeAddress();
  Int_t nmodules;
  fITS->InitModules(-1,nmodules); //-1->number of modules taken from AliITSgeom class kept in fITSgeom
                                  //nmodules is set

  fITS->FillModules(hits,0);

  for (Int_t imod=0; imod<240; ++imod){
    AliITSmodule *module = fITS->GetModule(imod);
    TObjArray *arrHits = module->GetHits();
    Int_t nhits = arrHits->GetEntriesFast();
    if (imod<80) {
      fAliITSQADataMakerSim->GetHitsData(fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(0.5,nhits);
      fAliITSQADataMakerSim->GetHitsData(1+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(imod,nhits);
    } else {
      fAliITSQADataMakerSim->GetHitsData(fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(1,nhits);
      fAliITSQADataMakerSim->GetHitsData(2+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(imod,nhits);
    }
    for (Int_t iHit=0; iHit<nhits; ++iHit) {
      AliITShit *hit = (AliITShit*) arrHits->At(iHit);
      Double_t xl,yl,zl,xl0,yl0,zl0;
      Double_t tof,tof0;
      hit->GetPositionL(xl,yl,zl,tof);
      hit->GetPositionL0(xl0,yl0,zl0,tof0);
      Float_t dyloc=TMath::Abs(yl-yl0)*10000.;
      fAliITSQADataMakerSim->GetHitsData(3+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(dyloc);
      Float_t edep=hit->GetIonization()*1000000;
      if(dyloc>180.){
        fAliITSQADataMakerSim->GetHitsData(4+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(edep);
      }
    }
  }
  return rv ; 
}


//_______________________________________________________________

Int_t AliITSQASPDDataMakerSim::GetOffset(AliQAv1::TASKINDEX_t task){
  // Returns histogram offset according to the specified task
  Int_t offset=0;
  if( task == AliQAv1::kHITS){
    offset=fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()];
  }
  else if( task == AliQAv1::kSDIGITS) {
    offset=fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()];
  }
  else if( task == AliQAv1::kDIGITS) {
    offset=fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()];
  }
  else {
    AliInfo("No task has been selected. TaskHisto set to zero.\n");
  }

  return offset;
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerSim::SetOffset(AliQAv1::TASKINDEX_t task, Int_t offset,Int_t specie ){
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

//_______________________________________________________________

Int_t AliITSQASPDDataMakerSim::GetTaskHisto(AliQAv1::TASKINDEX_t task) {
  // Returns the number of booked histograms for the selected task
  Int_t histotot=0;
  if( task == AliQAv1::kHITS) {
    histotot=fSPDhHTask ;
  }
  else if( task == AliQAv1::kSDIGITS) {
    histotot=fSPDhSTask;
  }
  else if( task == AliQAv1::kDIGITS) {
    histotot=fSPDhDTask ;
  }
  else {
    AliInfo("No task has been selected. TaskHisto set to zero.\n");
  }
  return histotot;

}
