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


/* $Id$ */

//---
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.
//  A. Mastroserio
//---

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TDirectory.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliT0digit.h" 
#include "AliT0hit.h"
#include "AliT0RecPoint.h"
#include "AliT0QADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliT0RawReader.h"

ClassImp(AliT0QADataMakerRec)
           
//____________________________________________________________________________ 
  AliT0QADataMakerRec::AliT0QADataMakerRec() : 
  AliQADataMakerRec(AliQA::GetDetName(AliQA::kT0), "T0 Quality Assurance Data Maker")

{
  // ctor
  /*
  for(Int_t i=0; i<24; i++) {
    fhHitsTime[i]=0x0;
   fhDigCFD[i]=0x0;
    fhDigLEDamp[i]=0x0;
    fhRecCFD[i]=0x0;
    fhRecLEDamp[i]=0x0;
    fhRecQTC[i]=0x0;
  }
  */
 //   fDetectorDir = fOutput->GetDirectory(GetName()) ;  
//   if (!fDetectorDir) 
//     fDetectorDir = fOutput->mkdir(GetName()) ;  
}

//____________________________________________________________________________ 
AliT0QADataMakerRec::AliT0QADataMakerRec(const AliT0QADataMakerRec& qadm) :
  AliQADataMakerRec() 
{
  //copy ctor 
  /*
  for(Int_t i=0; i<24; i++) {
    fhHitsTime[i]=0x0;
    fhDigCFD[i]=0x0;
    fhDigLEDamp[i]=0x0;
    fhRecCFD[i]=0x0;
    fhRecLEDamp[i]=0x0;
    fhRecQTC[i]=0x0;
  }
  */
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliT0QADataMakerRec& AliT0QADataMakerRec::operator = (const AliT0QADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliT0QADataMakerRec();
  new(this) AliT0QADataMakerRec(qadm);
  return *this;
}
//____________________________________________________________________________
void AliT0QADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kT0, task, list) ;
}

//____________________________________________________________________________
void AliT0QADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle

}
 
//____________________________________________________________________________ 
void AliT0QADataMakerRec::InitRaws()
{
  // create Raw histograms in Raw subdir
  printf("   AliT0QADataMakerRec::InitRaws() started\n");
  TString timename, ampname, qtcname;

  TH1F *fhRawCFD[24]; TH1F * fhRawLEDamp[24]; TH1F *fhRawQTC[24];

  for (Int_t i=0; i<24; i++)
    {
      timename ="hRawCFD";
      ampname = "hRawLED";
      qtcname = "hRawQTC";
      timename += i;
      ampname += i;
      qtcname += i;
      fhRawCFD[i] = new TH1F(timename.Data(), timename.Data(),100,100,5000);
      Add2RawsList( fhRawCFD[i],i);
      fhRawLEDamp[i] = new TH1F(ampname.Data(), ampname.Data(),100,120000,150000);
      Add2RawsList( fhRawLEDamp[i],i+24);
      fhRawQTC[i] = new TH1F(qtcname.Data(), qtcname.Data(),100,100,500);
      Add2RawsList( fhRawQTC[i],i+48);
     }
  
  TH1F* fhRawMean = new TH1F("hRawMean","online mean signal", 100,500,600);
  Add2RawsList( fhRawMean,72);
  
}

//____________________________________________________________________________ 

void AliT0QADataMakerRec::InitRecPoints()
{
  // create cluster histograms in RecPoint subdir
  /* 
 TH2F * fhRecCFD = new TH2F("fhRecCFD", " CFD reconstructed",25,-0.5,24.5,100,12,13);
  Add2DigitsList( fhRecCFD,0);
  TH2F *fhRecLEDamp = new TH2F("fhRecLEDamp", " amplitude LED reconstructed",25,-0.5,24.5,100,1000,1000);
  Add2DigitsList( fhRecLEDamp,1);
 TH2F * fhRecQTC = new TH2F("fhRecQTC", " amplitude QTC reconstructed",25,-0.5,24.5,100,1000,1000);
  Add2DigitsList( fhRecQTC,2);
 TH1F * fhRecMean = new TH1F("hRecMean"," reconstructed mean signal",100,500,600);
  Add2DigitsList( fhRecMean,3);
  */ 
  
  TString timename,ampname, qtcname;
  TH1F *fhRecCFD[24]; TH1F *fhRecLEDAmp[24];  TH1F * fhRecQTC[24];
  for (Int_t i=0; i<24; i++)
    {
      timename ="hRecCFD";
      ampname = "hRecLED";
      qtcname = "hRecQTC";
      timename += i;
      ampname += i;
      qtcname += i;
      fhRecCFD[i] = new TH1F(timename.Data(), timename.Data(),100,0,1000);
     Add2RecPointsList ( fhRecCFD[i],i);
      fhRecLEDAmp[i] = new TH1F(ampname.Data(), ampname.Data(),100,0,200);
    Add2RecPointsList ( fhRecLEDAmp[i],i+24);
      fhRecQTC[i] = new TH1F(qtcname.Data(), qtcname.Data(),100,0,200);
    Add2RecPointsList ( fhRecQTC[i],i+48);
     }
   
  TH1F *fhRecEff = new TH1F("hRecEff","Efficiency rec.points",25,-0.5,24.5);
 Add2RecPointsList ( fhRecEff,72);
  TH1F * fhRecMean = new TH1F("hRecMean"," reconstructed mean signal",100,500,600);
 Add2RecPointsList( fhRecMean,73);
  
}

//____________________________________________________________________________
void AliT0QADataMakerRec::InitESDs()
{
  //create ESDs histograms in ESDs subdir
  TH1F *fhESDMean = new TH1F("hESDmean"," ESD mean",100,0,100);
  Add2ESDsList(fhESDMean, 0) ;
 TH1F * fhESDVertex = new TH1F("hESDvertex","EAD vertex",100,-50,50);
  Add2ESDsList(fhESDVertex, 1) ;


}

//____________________________________________________________________________
void AliT0QADataMakerRec::MakeRaws( AliRawReader* rawReader)
{
  Int_t allData[110][5];
  for (Int_t i0=0; i0<105; i0++)
    {
      for (Int_t j0=0; j0<5; j0++) allData[i0][j0]=0;
    }
  //fills QA histos for RAW
AliT0RawReader *start = new AliT0RawReader(rawReader);
 
 while (rawReader->NextEvent()) {
   start->Next();
   for (Int_t i=0; i<105; i++) 
     for (Int_t iHit=0; iHit<5; iHit++)
       allData[i][iHit]= start->GetData(i,iHit);
   
   
   for (Int_t ik = 0; ik<24; ik+=2){
     for (Int_t iHt=0; iHt<5; iHt++){
       Int_t cc = ik/2;
       if(allData[cc+1][iHt]!=0){
	 GetRawsData(cc) -> Fill(allData[cc+1][iHt]-allData[0][0]);
	if(allData[ik+25][iHt]!=0 && allData[ik+26][iHt]!=0)
	  GetRawsData(cc+48)->Fill(allData[ik+26][iHt]-allData[ik+25][iHt]);
	if(allData[cc+13][iHt]!=0 )
	  GetRawsData(cc+24)->Fill(allData[cc+13][iHt]-allData[cc+1][iHt]);
	}
     }
   }
    
   for (Int_t ik = 24; ik<48; ik+=2) {
     for (Int_t iHt=0; iHt<5; iHt++) {
       Int_t cc = ik/2;
       if(allData[cc+45][iHt]!=0) {
	 GetRawsData(cc)->Fill(allData[cc+1][iHt]-allData[0][0]);
	if(allData[ik+57][iHt]!=0 && allData[ik+58][iHt]!=0)
	  GetRawsData(cc+48)->Fill(allData[ik+57][iHt]-allData[ik+58][iHt]);
	if(allData[cc+57][iHt]!=0 )
	  GetRawsData(cc+48)->Fill(allData[cc+57][iHt]-allData[cc+45][iHt]);
	}
     }
   }
   
 }

}

//____________________________________________________________________________
void AliT0QADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  //fills QA histos for clusters

   AliT0RecPoint* frecpoints= new AliT0RecPoint ();
    if (!frecpoints) {
    AliError("Reconstruct Fill ESD >> no recpoints found");
    return;
  }
  TBranch *brRec =clustersTree ->GetBranch("T0");
  if (brRec) {
    brRec->SetAddress(&frecpoints);
  }else{
      AliError(Form("EXEC Branch T0 rec not found "));
      return;
  } 
    
  brRec->GetEntry(0);
  
  for ( Int_t i=0; i<24; i++) {
    GetRecPointsData(i) -> Fill(frecpoints -> GetTime(i)); 
    GetRecPointsData(i+24) -> Fill(frecpoints -> GetAmp(i));
    GetRecPointsData(i+48) -> Fill(frecpoints->AmpLED(i));
    //  if(frecpoints -> GetTime(i) > 0) fhRecEff->Fill(i);
  }
     GetRecPointsData(72) ->Fill(frecpoints->GetMeanTime());
  
}

//____________________________________________________________________________
void AliT0QADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //fills QA histos for ESD

  GetESDsData(0) -> Fill(esd->GetT0());
  GetESDsData(1)-> Fill(esd->GetT0zVertex());

}

