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
#include "AliQA.h"
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
fSPDhDigits(0),
fSPDhSDigits(0),
fSPDhHits(0),
fGenOffset(0)
{
  //ctor used to discriminate OnLine-Offline analysis   
}

//____________________________________________________________________________ 
AliITSQASPDDataMakerSim::AliITSQASPDDataMakerSim(const AliITSQASPDDataMakerSim& qadm) :
TObject(),
fAliITSQADataMakerSim(qadm.fAliITSQADataMakerSim),
fSPDhDigits(qadm.fSPDhDigits),
fSPDhSDigits(qadm.fSPDhSDigits),
fSPDhHits(qadm.fSPDhHits),
fGenOffset(qadm.fGenOffset)
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
  AliDebug(1,"AliITSQADM::Start of SPD Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX_t /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  
  //AliQAChecker::Instance()->Run( AliQA::kITS , task, list);
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerSim::InitDigits()
{ 
  // Initialization for DIGIT data - SPD -
  fGenOffset = (fAliITSQADataMakerSim->fDigitsQAList)->GetEntries();
  //fSPDhDigits must be incremented by one unit every time a histogram is ADDED to the QA List

  Char_t name[50];
  Char_t title[50];

  TH1F *hlayer = new TH1F("LayPattern_SPD","Layer map - SPD",6,0.,6.);
  hlayer->GetXaxis()->SetTitle("Layer number");
  hlayer->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerSim->Add2DigitsList(hlayer,fGenOffset);
  fSPDhDigits++;
  
  TH1F **hmod = new TH1F*[2];
  for (Int_t iLay=0; iLay<2; iLay++) {
    sprintf(name,"ModPattern_SPD%d",iLay+1);
    sprintf(title,"Module map - SPD Layer %d",iLay+1);
    hmod[iLay]=new TH1F(name,title,240,0,240);
    hmod[iLay]->GetXaxis()->SetTitle("Module number");
    hmod[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerSim->Add2DigitsList(hmod[iLay],1+iLay+fGenOffset);
    fSPDhDigits++;
  }
  
  TH1F *hcolumns = new TH1F("Columns_SPD","Columns - SPD",160,0.,160.);
  hcolumns->GetXaxis()->SetTitle("Column number");
  hcolumns->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerSim->Add2DigitsList(hcolumns,3+fGenOffset);
  fSPDhDigits++;

  TH1F *hrows = new TH1F("Rows_SPD","Rows - SPD",256,0.,256.);
  hrows->GetXaxis()->SetTitle("Row number");
  hrows->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerSim->Add2DigitsList(hrows,4+fGenOffset);
  fSPDhDigits++;

  TH1F** hMultSPDdigits = new TH1F*[2];
  for (Int_t iLay=0; iLay<2; ++iLay) {
    sprintf(name,"DigitMultiplicity_SPD%d",iLay+1);
    sprintf(title,"Digit multiplicity - SPD Layer %d",iLay+1);
    hMultSPDdigits[iLay]=new TH1F(name,title,200,0.,200.);
    hMultSPDdigits[iLay]->GetXaxis()->SetTitle("Digit multiplicity");
    hMultSPDdigits[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerSim->Add2DigitsList(hMultSPDdigits[iLay], 5+iLay+fGenOffset);
    fSPDhDigits++;
  }

  TH2F *hMultSPDdig2MultSPDdig1 = new TH2F("DigitMultCorrelation_SPD","Digit multiplicity correlation - SPD",200,0.,200.,200,0.,200.);
  hMultSPDdig2MultSPDdig1->GetXaxis()->SetTitle("Digit multiplicity (Layer 1)");
  hMultSPDdig2MultSPDdig1->GetYaxis()->SetTitle("Digit multiplicity (Layer 2)");
  fAliITSQADataMakerSim->Add2DigitsList(hMultSPDdig2MultSPDdig1,7+fGenOffset);
  fSPDhDigits++;

  AliDebug(1,Form("%d SPD Digits histograms booked\n",fSPDhDigits));

}

//____________________________________________________________________________
void AliITSQASPDDataMakerSim::MakeDigits(TTree *digits)
{ 
  // Fill QA for DIGIT - SPD -
  AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  fITS->SetTreeAddress();
  TClonesArray *iITSdigits  = fITS->DigitsAddress(0);  // 0->SPD

  Int_t nDigitsL1=0;
  Int_t nDigitsL2=0;

  for (Int_t imod=0; imod<240; ++imod){
    digits->GetEvent(imod);
    Int_t ndigits = iITSdigits->GetEntries();
    if (imod<80) {
      fAliITSQADataMakerSim->GetDigitsData(0+fGenOffset)->Fill(0.5,ndigits);
      fAliITSQADataMakerSim->GetDigitsData(1+fGenOffset)->Fill(imod,ndigits);
      nDigitsL1+=ndigits;
    }
    else {
      fAliITSQADataMakerSim->GetDigitsData(0+fGenOffset)->Fill(1,ndigits);
      fAliITSQADataMakerSim->GetDigitsData(2+fGenOffset)->Fill(imod,ndigits);
      nDigitsL2+=ndigits;
    }
    for (Int_t idig=0; idig<ndigits; ++idig) {
      AliITSdigit *dig=(AliITSdigit*)iITSdigits->UncheckedAt(idig);
      Int_t col=dig->GetCoord1();  // cell number z
      Int_t row=dig->GetCoord2();  // cell number x
      fAliITSQADataMakerSim->GetDigitsData(3+fGenOffset)->Fill(col);
      fAliITSQADataMakerSim->GetDigitsData(4+fGenOffset)->Fill(row);
    }
  }
  fAliITSQADataMakerSim->GetDigitsData(5+fGenOffset)->Fill(nDigitsL1);
  fAliITSQADataMakerSim->GetDigitsData(6+fGenOffset)->Fill(nDigitsL2);
  fAliITSQADataMakerSim->GetDigitsData(7+fGenOffset)->Fill(nDigitsL1,nDigitsL2);
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerSim::InitSDigits()
{ 
  // Initialization for SDIGIT data - SPD -
  fGenOffset = (fAliITSQADataMakerSim->fSDigitsQAList)->GetEntries();
  //printf("--W-- AliITSQASPDDataMakerSim::InitSDigits()  fGenOffset= %d \n",fGenOffset);
  //fSPDhSDigits must be incremented by one unit every time a histogram is ADDED to the QA List
  
  Char_t name[50];
  Char_t title[50];

  TH1F *hlayer = new TH1F("LayPattern_SPD","Layer map - SPD",6,0.,6.);
  hlayer->GetXaxis()->SetTitle("Layer number");
  hlayer->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerSim->Add2SDigitsList(hlayer,fGenOffset);
  fSPDhSDigits++;

  TH1F **hmod = new TH1F*[2];
  for (Int_t iLay=0; iLay<2; ++iLay) {
    sprintf(name,"ModPattern_SPD%d",iLay+1);
    sprintf(title,"Module map - SPD Layer %d",iLay+1);
    hmod[iLay]=new TH1F(name,title,240,0,240);
    hmod[iLay]->GetXaxis()->SetTitle("Module number");
    hmod[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerSim->Add2SDigitsList(hmod[iLay],1+iLay+fGenOffset);
    fSPDhSDigits++;
  }
  

  AliDebug(1,Form("%d SPD SDigits histograms booked\n",fSPDhSDigits));

}

//____________________________________________________________________________
void AliITSQASPDDataMakerSim::MakeSDigits(TTree *sdigits)
{ 
  // Fill QA for SDIGIT - SPD -
  TBranch *brchSDigits = sdigits->GetBranch("ITS");
  for (Int_t imod=0; imod<240; ++imod){
    TClonesArray * sdig = new TClonesArray( "AliITSpListItem",1000 );
    brchSDigits->SetAddress( &sdig );
    brchSDigits->GetEvent(imod);
    Int_t nsdig=sdig->GetEntries();
    if (imod<80) {
      fAliITSQADataMakerSim->GetSDigitsData(0+fGenOffset)->Fill(0.5,nsdig);
      fAliITSQADataMakerSim->GetSDigitsData(1+fGenOffset)->Fill(imod,nsdig);
    }
    else {
      fAliITSQADataMakerSim->GetSDigitsData(0+fGenOffset)->Fill(1,nsdig);
      fAliITSQADataMakerSim->GetSDigitsData(2+fGenOffset)->Fill(imod,nsdig);
    }
    delete sdig;
  }

}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerSim::InitHits()
{ 
  // Initialization for HITS data - SPD -
  fGenOffset = (fAliITSQADataMakerSim->fHitsQAList)->GetEntries();
  //printf("--W-- AliITSQASPDDataMakerSim::InitHits()  fGenOffset= %d \n",fGenOffset);
  //fSPDhHits must be incremented by one unit every time a histogram is ADDED to the QA List
  Char_t name[50];
  Char_t title[50];

  TH1F *hlayer = new TH1F("LayPattern_SPD","Layer map - SPD",6,0.,6.);
  hlayer->GetXaxis()->SetTitle("Layer number");
  hlayer->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerSim->Add2HitsList(hlayer,fGenOffset);
  fSPDhHits++;

  TH1F **hmod = new TH1F*[2];
  for (Int_t iLay=0; iLay<2; ++iLay) {
    sprintf(name,"ModPattern_SPD%d",iLay+1);
    sprintf(title,"Module map - SPD Layer %d",iLay+1);
    hmod[iLay]=new TH1F(name,title,240,0,240);
    hmod[iLay]->GetXaxis()->SetTitle("Module number");
    hmod[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerSim->Add2HitsList(hmod[iLay],1+iLay+fGenOffset);
    fSPDhHits++;
  }

  TH1F *hhitlenght = new TH1F("Lenght_SPD","Hit lenght along y_{loc} coord",210,0.,210.);
  hhitlenght->GetXaxis()->SetTitle("Hit lenght [#mum]");
  hhitlenght->GetYaxis()->SetTitle("# hits");
  fAliITSQADataMakerSim->Add2HitsList(hhitlenght,3+fGenOffset);
  fSPDhHits++;

  TH1F *hEdepos = new TH1F("EnergyDeposit_SPD","Deposited energy distribution (y_{loc}>180 #mum)",150,0.,300.); 
  hEdepos->GetXaxis()->SetTitle("Deposited energy [keV]"); 
  hEdepos->GetYaxis()->SetTitle("# hits");
  fAliITSQADataMakerSim->Add2HitsList(hEdepos,4+fGenOffset);
  fSPDhHits++;

  AliDebug(1,Form("%d SPD Hits histograms booked\n",fSPDhHits));

}

//____________________________________________________________________________
void AliITSQASPDDataMakerSim::MakeHits(TTree *hits)
{ 
  // Fill QA for HITS - SPD -
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
      fAliITSQADataMakerSim->GetHitsData(fGenOffset)->Fill(0.5,nhits);
      fAliITSQADataMakerSim->GetHitsData(1+fGenOffset)->Fill(imod,nhits);
    } else {
      fAliITSQADataMakerSim->GetHitsData(fGenOffset)->Fill(1,nhits);
      fAliITSQADataMakerSim->GetHitsData(2+fGenOffset)->Fill(imod,nhits);
    }
    for (Int_t iHit=0; iHit<nhits; ++iHit) {
      AliITShit *hit = (AliITShit*) arrHits->At(iHit);
      Double_t xl,yl,zl,xl0,yl0,zl0;
      Double_t tof,tof0;
      hit->GetPositionL(xl,yl,zl,tof);
      hit->GetPositionL0(xl0,yl0,zl0,tof0);
      Float_t dyloc=TMath::Abs(yl-yl0)*10000.;
      fAliITSQADataMakerSim->GetHitsData(3+fGenOffset)->Fill(dyloc);
      Float_t edep=hit->GetIonization()*1000000;
      if(dyloc>180.){
        fAliITSQADataMakerSim->GetHitsData(4+fGenOffset)->Fill(edep);
      }
    }
  }
}
