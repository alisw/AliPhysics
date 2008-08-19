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
void AliT0QADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list)
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
  //  printf("   AliT0QADataMakerRec::InitRaws() started\n");
  TString timename, ampname, qtcname;
  TString timeCalname, ampCalname, qtcCalname;

  TH1F* fhRefPoint = new TH1F("hRefPoint","Ref Point", 10,1252170, 1252180);
  Add2RawsList( fhRefPoint,0);
   
  TH1F *fhRawCFD[24]; TH1F * fhRawLEDamp[24]; TH1F *fhRawQTC[24];
  TH1F *fhRawCFDcal[24]; TH1F * fhRawLEDampcal[24]; TH1F *fhRawQTCcal[24];
  
  for (Int_t i=0; i<24; i++)
    {
      timename ="hRawCFD";
      ampname = "hRawLED";
      qtcname = "hRawQTC";
      timename += i;
      ampname += i;
      qtcname += i;
      fhRawCFD[i] = new TH1F(timename.Data(), timename.Data(),500,2100,2800);
      Add2RawsList( fhRawCFD[i],i+1);
      fhRawLEDamp[i] = new TH1F(ampname.Data(), ampname.Data(),100,300,600);
      Add2RawsList( fhRawLEDamp[i],i+24+1);
      fhRawQTC[i] = new TH1F(qtcname.Data(), qtcname.Data(),1500,1000,7000);
      Add2RawsList( fhRawQTC[i],i+48+1);
     }
  const Bool_t saveForCorr = kTRUE;
  TH1F* fhRawMean = new TH1F("hRawMean","online mean signal", 100,2400,2500);
  Add2RawsList( fhRawMean,73, saveForCorr);
  TH1F* fhRawVertex = new TH1F("hRawVertex","online vertex signal", 100,0,600);
  Add2RawsList( fhRawVertex,74,saveForCorr );
  TH1F* fhRawORA = new TH1F("hRawORA","online OR A", 100,2500,2800);
  Add2RawsList( fhRawORA,75, saveForCorr);
  TH1F* fhRawORC = new TH1F("hRawORC","online OR C", 100,2000,2300);
  Add2RawsList( fhRawORC,76);
  for (Int_t i=0; i<24; i++)
    {
      // for events with trigger CALIBRATION_EVENT
      timeCalname ="hRawCFDcal";
      ampCalname = "hRawLEDcal";
      qtcCalname = "hRawQTCcal";
      timeCalname += i;
      ampCalname += i;
      qtcCalname += i;
      fhRawCFDcal[i] = new TH1F(timeCalname.Data(), timeCalname.Data(),10000,0,10000);
      Add2RawsList( fhRawCFDcal[i],76+i+1);
      fhRawLEDampcal[i] = new TH1F(ampCalname.Data(), ampCalname.Data(),100,300,600);
      Add2RawsList( fhRawLEDampcal[i],76+i+24+1);
      fhRawQTCcal[i] = new TH1F(qtcCalname.Data(), qtcCalname.Data(),1000,0,7000);
      Add2RawsList( fhRawQTCcal[i],76+i+48+1);
    }
  TH1F* fhRawMeanCal = new TH1F("hRawMeanCal","online mean signal, calibration event",
				10000,0,10000);
  Add2RawsList( fhRawMeanCal,149, saveForCorr);
  TH1F* fhRawVertexCal = new TH1F("hRawVertexCal","online vertex signal, calibration even ",
				  10000,0,10000);
  Add2RawsList( fhRawVertexCal,150, saveForCorr);
  TH1F* fhRawORAcal = new TH1F("hRawORAcal","online OR A", 10000,0,10000);
  Add2RawsList( fhRawORAcal,151, saveForCorr );
  TH1F* fhRawORCcal = new TH1F("hRawORCcal","online OR C", 10000,0,10000);
  Add2RawsList( fhRawORCcal,152, saveForCorr);
  TH1F* fhMultcal = new TH1F("hMultcal","full mulltiplicity", 10000,0,10000);
  Add2RawsList( fhMultcal,153 );
  TH1F* fhMultScal = new TH1F("hMultScal","full multiplicity with semi-central trigger",
			      10000,0,10000);
  Add2RawsList( fhMultScal,154);
  TH1F* fhMultCcal = new TH1F("hMultCcal","full multiplicity with central trigger", 
			      1000,0,10000);
  Add2RawsList( fhMultCcal,155);
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
  
  //  printf(" !!!!!  AliT0QADataMakerRec::InitRecPoints() started\n");
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
      fhRecCFD[i] = new TH1F(timename.Data(), timename.Data(),100,2100,2800);
     Add2RecPointsList ( fhRecCFD[i],i);
      fhRecLEDAmp[i] = new TH1F(ampname.Data(), ampname.Data(),100,0, 100);
    Add2RecPointsList ( fhRecLEDAmp[i],i+24);
      fhRecQTC[i] = new TH1F(qtcname.Data(), qtcname.Data(),100,0,100);
    Add2RecPointsList ( fhRecQTC[i],i+48);
     }
   
  TH1F *fhOnlineMean = new TH1F("hOnlineMean","online mean",100,2400,2500);
  Add2RecPointsList ( fhOnlineMean,72);
  TH1F * fhRecMean = new TH1F("hRecMean"," reconstructed mean signal",100,2400,2500);
  Add2RecPointsList( fhRecMean,73);
  //  printf(" !!!!!!  AliT0QADataMakerRec::InitRecPoints() ended\n");
  
}

//____________________________________________________________________________
void AliT0QADataMakerRec::InitESDs()
{
  //create ESDs histograms in ESDs subdir

  TH1F *fhESDMean = new TH1F("hESDmean"," ESD mean",100,2400,2500);
  Add2ESDsList(fhESDMean, 0) ;
  TH1F * fhESDVertex = new TH1F("hESDvertex","EAD vertex",100,-50,50);
  Add2ESDsList(fhESDVertex, 1) ;
  

}

//____________________________________________________________________________
void AliT0QADataMakerRec::MakeRaws( AliRawReader* rawReader)
{
	rawReader->Reset() ; 
  //fills QA histos for RAW
  Int_t shift=0;
  AliT0RawReader *start = new AliT0RawReader(rawReader);
  //  start->Next();
  if (! start->Next())
    AliDebug(1,Form(" no raw data found!!"));
  else
    {  

      UInt_t type =rawReader->GetType();
      //     cout<<" !!!!! new event type = "<<type<<endl;
      Int_t allData[110][5];
      for (Int_t i0=0; i0<105; i0++)
	{
	  for (Int_t j0=0; j0<5; j0++) allData[i0][j0]=0;
	}
      for (Int_t i=0; i<105; i++) 
	for (Int_t iHit=0; iHit<5; iHit++)
	  allData[i][iHit]= start->GetData(i,iHit);
      
      GetRawsData(0) -> Fill( allData[0][0]);
      allData[0][0] = allData[0][0] - 7000; 
      if (type == 8) shift=76;
      if (type == 7) shift=0;
	    
      for (Int_t ik = 0; ik<12; ik++){
	for (Int_t iHt=0; iHt<5; iHt++){
	  if(allData[ik+1][iHt]>0){
	      GetRawsData(shift+ik+1) -> Fill(allData[ik+1][iHt]-allData[0][0]);
	      //	      cout<<" type "<<type<<" shift "<<shift<<" koef "<<ik<<" index "<<shift+ik<<
	      //   	" time "<<allData[ik+1][iHt]-allData[0][0]<<endl;
	      if(allData[2*ik+25][iHt] > 0 && allData[2*ik+26][iHt] > 0)
		GetRawsData(shift+ik+48+1)->Fill(allData[2*ik+25][iHt]-allData[2*ik+26][iHt]);
	      if(allData[ik+13][iHt]!=0 )
		GetRawsData(shift+ik+24+1)->Fill(allData[ik+13][iHt]-allData[ik+1][iHt]);
	  
	}
      }
      }
      for (Int_t ik = 12; ik<24; ik++) {
	for (Int_t iHt=0; iHt<5; iHt++) {
	  if(allData[ik+45][iHt]>0) {
	    GetRawsData(shift+ik+1)->Fill(allData[ik+45][iHt]-allData[0][0]);
	    //  cout<<" type "<<type<<" shift "<<shift<<" index "<<shift+ik<<
	    // " time "<<allData[ik+1][iHt]-allData[0][0]<<endl;

	    if(allData[2*ik+57][iHt]!=0 && allData[2*ik+58][iHt]!=0)
	      GetRawsData(shift+ik+48+1)->Fill(allData[2*ik+57][iHt]-allData[2*ik+58][iHt]);
	    if(allData[ik+57][iHt] > 0 )
	      GetRawsData(shift+ik+24+1)->Fill(allData[ik+57][iHt]-allData[ik+45][iHt]);
	  }
	}
      }
      

      if(type == 7)
	{
	  for (Int_t iHt=0; iHt<5; iHt++) {
	    GetRawsData(73)->Fill(allData[49][iHt]-allData[0][0]);
	    GetRawsData(74)->Fill(allData[50][iHt]-allData[0][0]);
	    GetRawsData(75)->Fill(allData[51][iHt]-allData[0][0]);
	    GetRawsData(76)->Fill(allData[52][iHt]-allData[0][0]);
	  }
	}
      if(type == 8)
	{
	  for (Int_t iHt=0; iHt<5; iHt++) {
	    GetRawsData(149)->Fill(allData[49][iHt]-allData[0][0]);
	    GetRawsData(150)->Fill(allData[50][iHt]-allData[0][0]);
	    GetRawsData(151)->Fill(allData[51][iHt]-allData[0][0]);
	    GetRawsData(152)->Fill(allData[52][iHt]-allData[0][0]);
	    /*	  cout<<" and "<<allData[49][0]-allData[0][0]<<
	    " vertex "<<allData[50][0]-allData[0][0]<<
	    " ORA "<<allData[51][0]-allData[0][0]<<
	    " ORC "<<allData[52][0]-allData[0][0]<<endl;*/
	    
	    GetRawsData(153)->Fill(allData[53][iHt]-allData[54][iHt]);
	    if(allData[55][iHt])  GetRawsData(154)->Fill(allData[53][iHt]-allData[54][iHt]);
	    if(allData[55][iHt])  GetRawsData(155)->Fill(allData[53][iHt]-allData[54][iHt]);
	  }
	}  
      delete start;
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
     GetRecPointsData(72) ->Fill(frecpoints->GetOnlineMean());
     GetRecPointsData(73) ->Fill(frecpoints->GetMeanTime());
  
}

//____________________________________________________________________________
void AliT0QADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //fills QA histos for ESD

  GetESDsData(0) -> Fill(esd->GetT0());
  GetESDsData(1)-> Fill(esd->GetT0zVertex());

}

