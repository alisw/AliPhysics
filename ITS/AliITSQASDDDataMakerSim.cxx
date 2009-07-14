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

/* $Id$ */

//  *************************************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  contained in a DB
//  -------------------------------------------------------------
//  W. Ferrarese + P. Cerello Feb 2008
//  INFN Torino

// --- ROOT system ---
#include <TTree.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQASDDDataMakerSim.h"
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliQADataMakerSim.h"
#include "AliITSQADataMakerSim.h"
#include "AliRawReader.h"
#include "AliITSdigit.h"
#include "AliITS.h"
#include "AliITSmodule.h"
#include "AliITShit.h"
#include "AliITSLoader.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSpList.h"

ClassImp(AliITSQASDDDataMakerSim)

//____________________________________________________________________________ 
AliITSQASDDDataMakerSim::AliITSQASDDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim) :
TObject(),
fAliITSQADataMakerSim(aliITSQADataMakerSim),
fSDDhHTask(0),
fSDDhSTask(0),
fSDDhDTask(0),
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
AliITSQASDDDataMakerSim::AliITSQASDDDataMakerSim(const AliITSQASDDDataMakerSim& qadm) :
TObject(),
fAliITSQADataMakerSim(qadm.fAliITSQADataMakerSim),
fSDDhHTask(qadm.fSDDhHTask),
fSDDhSTask(qadm.fSDDhSTask),
fSDDhDTask(qadm.fSDDhDTask),
fGenOffsetH(qadm.fGenOffsetH),
fGenOffsetS(qadm.fGenOffsetS),
fGenOffsetD(qadm.fGenOffsetD)
{
  //copy ctor 
  fAliITSQADataMakerSim->SetName((const char*)qadm.fAliITSQADataMakerSim->GetName()) ; 
  fAliITSQADataMakerSim->SetTitle((const char*)qadm.fAliITSQADataMakerSim->GetTitle());

}

//__________________________________________________________________
AliITSQASDDDataMakerSim& AliITSQASDDDataMakerSim::operator = (const AliITSQASDDDataMakerSim& qac )
{
  // Equal operator.
  this->~AliITSQASDDDataMakerSim();
  new(this) AliITSQASDDDataMakerSim(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Start of SDD Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSDM instantiates checker with Run(AliQAv1::kITS, task, list)\n"); 
  //AliQAChecker::Instance()->Run( AliQAv1::kITS , task, list);
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerSim::InitDigits()
{ 
  // Initialization for DIGIT data - SDD -  
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ;
  Int_t rv = 0 ; 
  //fGenOffsetD = (fAliITSQADataMakerSim->fDigitsQAList[AliRecoParam::kDefault])->GetEntries();
  //fSDDhTask must be incremented by one unit every time a histogram is ADDED to the QA List
  TH1F* h0=new TH1F("SDD DIGITS Module Pattern","SDD DIGITS Module Pattern",260,239.5,499.5);       //hmod
  h0->GetXaxis()->SetTitle("SDD Module Number");
  h0->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerSim->Add2DigitsList(h0,fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSDDhDTask ++;
  TH1F* h1=new TH1F("SDD Anode Distribution","SDD DIGITS Anode Distribution",512,-0.5,511.5);      //hanocc
  h1->GetXaxis()->SetTitle("Anode Number");
  h1->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerSim->Add2DigitsList(h1,1+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSDDhDTask ++;
  TH1F* h2=new TH1F("SDD Tbin Distribution","SDD DIGITS Tbin Distribution",256,-0.5,255.5);      //htbocc
  h2->GetXaxis()->SetTitle("Tbin Number");
  h2->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerSim->Add2DigitsList(h2,2+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSDDhDTask ++;
  TH1F* h3=new TH1F("SDD ADC Counts Distribution","SDD DIGITS ADC Counts Distribution",200,0.,1024.);          //hsig
  h3->GetXaxis()->SetTitle("ADC Value");
  h3->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerSim->Add2DigitsList(h3,3+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSDDhDTask ++;
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Digits histograms booked\n",fSDDhDTask));
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASDDDataMakerSim::MakeDigits(TTree * digits)
{ 

  // Fill QA for DIGIT - SDD -
  Int_t rv = 0 ; 

  AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  fITS->SetTreeAddress();
  TClonesArray *iITSdigits  = fITS->DigitsAddress(1);
  for(Int_t i=0; i<260; i++){
    Int_t nmod=i+240;
    digits->GetEvent(nmod);
    Int_t ndigits = iITSdigits->GetEntries();
    fAliITSQADataMakerSim->GetDigitsData(fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(nmod,ndigits);
    for (Int_t idig=0; idig<ndigits; idig++) {
      AliITSdigit *dig=(AliITSdigit*)iITSdigits->UncheckedAt(idig);
      Int_t iz=dig->GetCoord1();  // cell number z
      Int_t ix=dig->GetCoord2();  // cell number x
      Int_t sig=dig->GetSignal();
      fAliITSQADataMakerSim->GetDigitsData(1+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(iz);
      fAliITSQADataMakerSim->GetDigitsData(2+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(ix);
      fAliITSQADataMakerSim->GetDigitsData(3+fGenOffsetD[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(sig);
    }
  }
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerSim::InitSDigits()
{ 
  // Initialization for SDIGIT data - SDD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ;
  Int_t rv = 0 ; 
  //fGenOffsetS = (fAliITSQADataMakerSim->fSDigitsQAList[AliRecoParam::kDefault])->GetEntries();
  //fSDDhTask must be incremented by one unit every time a histogram is ADDED to the QA List
  TH1F* h0=new TH1F("SDD SDIGITS Module Pattern","SDIGITS SDD Module Pattern",260,239.5,499.5);       //hmod
  h0->GetXaxis()->SetTitle("SDD Module Number");
  h0->GetYaxis()->SetTitle("# SDIGITS");
  rv = fAliITSQADataMakerSim->Add2SDigitsList(h0,fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSDDhSTask ++;
  TH1F* h1=new TH1F("SDD Anode Distribution","SDIGITS SDD Anode Distribution",512,-0.5,511.5);      //hanocc
  h1->GetXaxis()->SetTitle("Anode Number");
  h1->GetYaxis()->SetTitle("# SDIGITS");
  rv = fAliITSQADataMakerSim->Add2SDigitsList(h1,1+fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSDDhSTask ++;
  TH1F* h2=new TH1F("SDD Tbin Distribution","SDIGITS SDD Tbin Distribution",256,-0.5,255.5);      //htbocc
  h2->GetXaxis()->SetTitle("Tbin Number");
  h2->GetYaxis()->SetTitle("# SDIGITS");
  rv = fAliITSQADataMakerSim->Add2SDigitsList(h2,2+fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()]);
  fSDDhSTask ++;
  TH1F* h3=new TH1F("SDD ADC Counts Distribution","SDIGITS SDD ADC Counts Distribution",200,0.,1024.);          //hsig
  h3->GetXaxis()->SetTitle("ADC Value");
  h3->GetYaxis()->SetTitle("# SDIGITS");
  rv = fAliITSQADataMakerSim->Add2SDigitsList(h3,3+fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSDDhSTask ++;

  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD SDigits histograms booked\n",fSDDhSTask));
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASDDDataMakerSim::MakeSDigits(TTree * sdigits)
{ 
  // Fill QA for SDIGIT - SDD -
  Int_t rv = 0 ; 
  
   AliITSsegmentationSDD* seg = new AliITSsegmentationSDD();
  Int_t nan=seg->Npz();
  Int_t ntb=seg->Npx();
  Int_t scaleSize=4;
  AliITSpList* list=new AliITSpList(nan,ntb*scaleSize);

  //AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  //fITS->SetTreeAddress();
  //TClonesArray *ITSdigits  = fITS->DigitsAddress(1);
  //TFile *sper = new TFile("sper.root","CREATE"); //agginto a mano x prova
  //digits->Write();
  //sper->Close();


  TBranch *brchSDigits = sdigits->GetBranch("ITS");
  static TClonesArray * sdig ; 
  if (! sdig )
    sdig = new TClonesArray( "AliITSpListItem",1000 );
  for(Int_t id=0; id<260; id++){
    Int_t nmod=id+240;
    brchSDigits->SetAddress( &sdig );
    brchSDigits->GetEvent(nmod);
    Int_t nsdig=sdig->GetEntries();
    fAliITSQADataMakerSim->GetSDigitsData(fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(nmod,nsdig);
    for(Int_t i=0;i<nsdig;i++){
      AliITSpListItem *cell=(AliITSpListItem*)sdig->At(i);
      Float_t sig=cell->GetSignal();
      Int_t idx=cell->GetIndex();
      Int_t ia,it;
      list->GetCell(idx,ia,it);
      fAliITSQADataMakerSim->GetSDigitsData(1+fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(ia);
      fAliITSQADataMakerSim->GetSDigitsData(2+fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(it);
      fAliITSQADataMakerSim->GetSDigitsData(3+fGenOffsetS[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(sig);
    }
    sdig->Clear();
  }
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerSim::InitHits()
{ 

  // Initialization for HITS data - SDD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ;
  Int_t rv = 0 ; 

  //fGenOffsetH = (fAliITSQADataMakerSim->fHitsQAList[AliRecoParam::kDefault])->GetEntries();
  //fSDDhTask must be incremented by one unit every time a histogram is ADDED to the QA List
  //printf("AliITSQASDDDataMakerSim::InitHits called \n");
  TH1F *h0=new TH1F("SDD HITS Module Pattern","SDD HITS Module Pattern",260,239.5,499.5);  
  h0->GetXaxis()->SetTitle("SDD Module Number");
  h0->GetYaxis()->SetTitle("# HITS");
  rv = fAliITSQADataMakerSim->Add2HitsList(h0,fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSDDhHTask ++;
  TH1F *h1=new TH1F("SDD HIT lenght along local Y Coord","SDD HIT lenght along local Y Coord",200,0.,350.);
  h1->GetXaxis()->SetTitle("HIT lenght (um)");
  h1->GetYaxis()->SetTitle("# HITS");
  rv = fAliITSQADataMakerSim->Add2HitsList(h1,1+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSDDhHTask ++;
  TH1F *h2=new TH1F("SDD HIT lenght along local Y Coord - Zoom","SDD HIT lenght along local Y Coord - Zoom",200,250.,350.);
  h2->GetXaxis()->SetTitle("HIT lenght (um)");
  h2->GetYaxis()->SetTitle("# HITS");
  rv = fAliITSQADataMakerSim->Add2HitsList(h2,2+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSDDhHTask ++;
  TH1F *h3=new TH1F("SDD Deposited Energy Distribution (loc Y > 200um)","SDD HITS Deposited Energy Distribution (loc Y > 200um)",200,0.,350.);
  h3->GetXaxis()->SetTitle("ADC counts ");
  h3->GetYaxis()->SetTitle("# HITS");
  rv = fAliITSQADataMakerSim->Add2HitsList(h3,3+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()], !expert, image);
  fSDDhHTask ++;
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Hits histograms booked\n",fSDDhHTask));
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASDDDataMakerSim::MakeHits(TTree * hits)
{ 
  // Fill QA for HITS - SDD -
  Int_t rv = 0 ; 

   AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  fITS->SetTreeAddress();
  Int_t nmodules;
  if(!(fITS->InitModules(-1,nmodules))){
    AliError("ITS geometry not available - nothing done");
    return rv;
  }
 
  fITS->FillModules(hits,0);

  for(Int_t i=0; i<260; i++){
    Int_t nmod=i+240;
    AliITSmodule *modu = fITS->GetModule(nmod);
    TObjArray *arrHits = modu->GetHits();
    Int_t nhits = arrHits->GetEntriesFast();
    ////printf("--w--AliITSQASDDDataMakerSim::MakeHits  nhits = %d\n",nhits);
    for (Int_t iHit=0;iHit<nhits;iHit++) {
      AliITShit *hit = (AliITShit*) arrHits->At(iHit);
      fAliITSQADataMakerSim->GetHitsData(fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(nmod);
      Double_t xl,yl,zl,xl0,yl0,zl0;
      Double_t tof,tof0;
      hit->GetPositionL(xl,yl,zl,tof);
      hit->GetPositionL0(xl0,yl0,zl0,tof0);
      Float_t dyloc=TMath::Abs(yl-yl0)*10000.;
      fAliITSQADataMakerSim->GetHitsData(1+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(dyloc);
      Float_t edep=hit->GetIonization()*1000000;
      if(dyloc>200.){ 
        fAliITSQADataMakerSim->GetHitsData(2+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(edep);
        fAliITSQADataMakerSim->GetHitsData(3+fGenOffsetH[fAliITSQADataMakerSim->GetEventSpecie()])->Fill(dyloc);
      }
    }
  }
  return rv ; 
}

//_______________________________________________________________

Int_t AliITSQASDDDataMakerSim::GetOffset(AliQAv1::TASKINDEX_t task){
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
void AliITSQASDDDataMakerSim::SetOffset(AliQAv1::TASKINDEX_t task, Int_t offset,Int_t specie ){
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

Int_t AliITSQASDDDataMakerSim::GetTaskHisto(AliQAv1::TASKINDEX_t task) {
  // Returns the number of booked histograms for the selected task
  Int_t histotot=0;
  if( task == AliQAv1::kHITS) {
    histotot=fSDDhHTask ;  
  }
  else if( task == AliQAv1::kSDIGITS) {
    histotot=fSDDhSTask;   
  }
  else if( task == AliQAv1::kDIGITS) {
    histotot=fSDDhDTask ;   
  }
  else {
    AliInfo("No task has been selected. TaskHisto set to zero.\n");
  }
  return histotot;

}
