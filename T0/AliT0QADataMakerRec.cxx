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
//  Alla.Maevskaya@cern.ch
//---

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH2F.h> 
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

#include "Riostream.h"
ClassImp(AliT0QADataMakerRec)
           
//____________________________________________________________________________ 
  AliT0QADataMakerRec::AliT0QADataMakerRec() : 
AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kT0), 
		  "T0 Quality Assurance Data Maker"),
  fnEvent(0)

{
  // ctor
  for (Int_t i=0; i<6; i++) {
    fNumTriggers[i]=0;
    fNumTriggersCal[i]=0;
  }
  for (Int_t i=0; i<24; i++)
    {
      feffC[i]=0;
      feffA[i]=0;
    }
}

//____________________________________________________________________________ 
AliT0QADataMakerRec::AliT0QADataMakerRec(const AliT0QADataMakerRec& qadm) :
  AliQADataMakerRec(),
  fnEvent(0)
  
{
  //copy ctor 
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
void AliT0QADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kT0, task, list) ;
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    SetEventSpecie(AliRecoParam::ConvertIndex(specie)) ; 
    if ( task == AliQAv1::kRAWS ) {
      const Char_t *triggers[6] = {"mean", "vertex","ORA","ORC","central","semi-central"};
      for (Int_t itr=0; itr<6; itr++) {
        GetRawsData(197)->Fill(triggers[itr], fNumTriggersCal[itr]);
        GetRawsData(197)->SetBinContent(itr+1, fNumTriggersCal[itr]);
      }  
      GetRawsData(205)->SetOption("COLZ");
      GetRawsData(206)->SetOption("COLZ");
      GetRawsData(207)->SetOption("COLZ");
      GetRawsData(205)->GetXaxis()->SetTitle("#PMT");
      GetRawsData(206)->GetXaxis()->SetTitle("#PMT");
      GetRawsData(205)->GetYaxis()->SetTitle("NeventsReg/Nevents");
      GetRawsData(206)->GetYaxis()->SetTitle("NeventsReg/Nevents");
      GetRawsData(207)->GetXaxis()->SetTitle("#PMT");
      GetRawsData(207)->GetYaxis()->SetTitle("Charge, #channels");
    }
    if ( task == AliQAv1::kRECPOINTS) {
      GetRecPointsData(0)->SetOption("COLZ");
      GetRecPointsData(1)->SetOption("COLZ");
      GetRecPointsData(0)->GetXaxis()->SetTitle("#PMT");
      GetRecPointsData(1)->GetXaxis()->SetTitle("#PMT");
      GetRecPointsData(0)->GetYaxis()->SetTitle("CFD time");
      GetRecPointsData(1)->GetYaxis()->SetTitle("Charge, #channels");
    }
  }
}

//____________________________________________________________________________
void AliT0QADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  fnEvent=0;

}
 
//____________________________________________________________________________ 
void AliT0QADataMakerRec::InitRaws()
{
  // create Raw histograms in Raw subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TString timename, ampname, qtcname, ledname;
  TString timeCalname, ampCalname, ledCalname, qtcCalname;

  TH1F* fhRefPoint = new TH1F("hRefPoint","Ref Point", 10,1252170, 1252180);
  Add2RawsList( fhRefPoint,0, !expert, image, !saveCorr);
   
  TH1F *fhRawCFD[24]; TH1F * fhRawLEDamp[24];
  TH1F *fhRawQTC[24]; TH1F * fhRawLED[24];
  TH1F *fhRawCFDcal[24]; TH1F * fhRawLEDampcal[24]; 
  TH1F *fhRawQTCcal[24];  TH1F * fhRawLEDcal[24];
  
  for (Int_t i=0; i<24; i++)
    {
      timename ="hRawCFD";
      ledname = "hRawLED";
      qtcname = "hRawQTC";
      ampname = "hRawLEDminCFD";
      timename += i;
      ampname += i;
      qtcname += i;
      ledname += i;
      fhRawCFD[i] = new TH1F(timename.Data(), timename.Data(),500,-250,250);
      Add2RawsList( fhRawCFD[i],i+1, !expert, image, !saveCorr);
      fhRawLED[i] = new TH1F(ledname.Data(), ledname.Data(),2000,-1000,1000);
      Add2RawsList( fhRawLED[i],i+24+1, !expert, image, !saveCorr);
      fhRawLEDamp[i] = new TH1F(ampname.Data(), ampname.Data(),100,300,600);
      Add2RawsList( fhRawLEDamp[i],i+48+1, !expert, image, !saveCorr);
      fhRawQTC[i] = new TH1F(qtcname.Data(), qtcname.Data(),1500,1000,7000);
      Add2RawsList( fhRawQTC[i],i+72+1, !expert, image, !saveCorr);
     }
  TH1F* fhRawTrigger = new TH1F("hRawTrigger"," phys triggers",5,0,5);
  Add2RawsList(fhRawTrigger ,97, !expert, image, !saveCorr);
  
  TH1F* fhRawMean = new TH1F("hRawMean","online mean signal", 100,2400,2500);
  Add2RawsList( fhRawMean,98, !expert, image, !saveCorr);
  TH1F* fhRawVertex = new TH1F("hRawVertex","online vertex signal", 100,0,600);
  Add2RawsList( fhRawVertex,99, !expert, image, !saveCorr);
  TH1F* fhRawORA = new TH1F("hRawORA","online OR A", 100,2500,2800);
  Add2RawsList( fhRawORA,100, !expert, image, !saveCorr);
  TH1F* fhRawORC = new TH1F("hRawORC","online OR C", 100,2000,2300);
  Add2RawsList( fhRawORC,101, !expert, image, !saveCorr);
  
  for (Int_t i=0; i<24; i++)
    {
      // for events with trigger CALIBRATION_EVENT
      timeCalname ="hRawCFDcal";
      ledCalname = "hRawLEDcal";
      ampCalname = "hRawLEDminCFDcal";
      qtcCalname = "hRawQTCcal";
      timeCalname += i;
      ledCalname += i;
      ampCalname += i;
      qtcCalname += i;
      fhRawCFDcal[i] = new TH1F(timeCalname.Data(), timeCalname.Data(),2000,-1000,1000);
      Add2RawsList( fhRawCFDcal[i],101+i+1, !expert, image, !saveCorr);
      fhRawLEDcal[i] = new TH1F(ledCalname.Data(), ledCalname.Data(),2000,-1000,1000);
      Add2RawsList( fhRawLEDcal[i],101+i+24+1, !expert, image, !saveCorr);
      fhRawLEDampcal[i] = new TH1F(ampCalname.Data(), ampCalname.Data(),300,300,600);
      Add2RawsList( fhRawLEDampcal[i],101+i+48+1, !expert, image, !saveCorr);
      fhRawQTCcal[i] = new TH1F(qtcCalname.Data(), qtcCalname.Data(),1000,0,7000);
      Add2RawsList( fhRawQTCcal[i],101+i+72+1, !expert, image, !saveCorr);
    }

  TH1F* fhRawTriggerCal = new TH1F("hRawTriggerCal"," laser triggers",6,0,6);
  Add2RawsList(fhRawTriggerCal ,197 , !expert, !image, saveCorr);

  TH1F* fhRawMeanCal = new TH1F("hRawMeanCal","online mean signal, calibration event",
				10000,0,10000);
  Add2RawsList( fhRawMeanCal,198);
  TH1F* fhRawVertexCal = new TH1F("hRawVertexCal","online vertex signal, calibration even ",
				  10000,0,10000);
  Add2RawsList( fhRawVertexCal,199, !expert, image, !saveCorr);
  TH1F* fhRawORAcal = new TH1F("hRawORAcal","online OR A", 10000,0,10000);
  Add2RawsList( fhRawORAcal,200, !expert, image, !saveCorr );
  TH1F* fhRawORCcal = new TH1F("hRawORCcal","online OR C", 10000,0,10000);
  Add2RawsList( fhRawORCcal,201, !expert, image, !saveCorr);
  TH1F* fhMultcal = new TH1F("hMultcal","full mulltiplicity", 10000,0,10000);
  Add2RawsList( fhMultcal,202, !expert, image, !saveCorr );
  TH1F* fhMultScal = new TH1F("hMultScal","full multiplicity with semi-central trigger",
			      10000,0,10000);
  Add2RawsList( fhMultScal,203, !expert, image, !saveCorr);
  TH1F* fhMultCcal = new TH1F("hMultCcal","full multiplicity with central trigger", 
			      1000,0,10000);
  Add2RawsList( fhMultCcal,204, !expert, image, !saveCorr);

  //  TH2F* fhEffCFD = new TH2F("hEffCFD"," CFD time",24, 0 ,24, 
  //		      100,-500,500);
  TH2F* fhEffCFD = new TH2F("hEffCFD"," CFD time",24, 0 ,24, 50, 0,5); 

  Add2RawsList( fhEffCFD,205, !expert, !image, saveCorr);
  TH2F* fhEffLED = new TH2F("hEffLED","LED time",24, 0 ,24, 
			    100, 0, 5);	
  //100,-500,500);

  Add2RawsList( fhEffLED,206, !expert, !image, saveCorr);
  TH2F* fhEffQTC = new TH2F("hEffQTC","QTC amplitude",24, 0 ,24, 
			      100,0,7000);
  Add2RawsList( fhEffQTC,207, !expert, !image, saveCorr);
  //yeys guide line
  //  TH2F* fhLineQTC = new TH2F("hLineQTC","QTC amplitude boeder",24, 0 ,24, 
  //			      100,0,7000);
//  Add2RawsList( fhLineQTC,208, !expert, !image, saveCorr);
}

//____________________________________________________________________________ 
void AliT0QADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH2F * fhDigCFD = new TH2F("fhDigCFD", " CFD digits",25,-0.5,24.5,100,100,1000);
  Add2DigitsList( fhDigCFD,0);
  TH2F *fhDigLEDamp = new TH2F("fhDigLEDamp", " LED-CFD digits",25,-0.5,24.5,100,100,1000);
  Add2DigitsList( fhDigLEDamp,1, !expert, image);
  TH2F * fhDigQTC = new TH2F("fhDigQTC", " QTC digits",25,-0.5,24.5,100,100,1000);
  Add2DigitsList( fhDigQTC,2, !expert, image);
  
  
  
}

//____________________________________________________________________________ 

void AliT0QADataMakerRec::InitRecPoints()
{
  // create cluster histograms in RecPoint subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH2F* fhRecCFD = new TH2F("hRecCFD"," CFD time",24, 0 ,24, 
			      100,-50,50);
  Add2RecPointsList ( fhRecCFD,0, !expert, image);

  TH2F* fhRecAmpDiff = new TH2F("hRecAmpDiff"," LED-CFD  min QTC amplitude",
				24, 0 ,24, 200,-10,10);
  Add2RecPointsList (fhRecAmpDiff, 1, !expert, image);
  
  TH1F *fhMean = new TH1F("hMean","online - rec mean",1000, -5000, 5000);
  Add2RecPointsList ( fhMean,2, !expert, image);
 }

//____________________________________________________________________________
void AliT0QADataMakerRec::InitESDs()
{
  //create ESDs histograms in ESDs subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *fhESDMean = new TH1F("hESDmean"," ESD mean",100,2400,2500);
  Add2ESDsList(fhESDMean, 0, !expert, image) ;
  TH1F * fhESDVertex = new TH1F("hESDvertex","ESD vertex",82,-30,30);
  Add2ESDsList(fhESDVertex, 1, !expert, image) ;
  

}

//____________________________________________________________________________
void AliT0QADataMakerRec::MakeRaws( AliRawReader* rawReader)
{
	rawReader->Reset() ; 
  //fills QA histos for RAW
  Int_t shift=0;
  Float_t effic ;

  AliT0RawReader *start = new AliT0RawReader(rawReader);
  //  start->Next();
  if (! start->Next())
    AliDebug(AliQAv1::GetQADebugLevel(),Form(" no raw data found!!"));
  else
    {  
      fnEvent++;

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
      
      if (allData[0][0]>0)    GetRawsData(0) -> Fill( allData[0][0]);
      //  allData[0][0] = allData[0][0] - 7000; 
      if (type == 8) shift=101;
      if (type == 7) shift=0;
      
      for (Int_t ik = 0; ik<12; ik++){
	for (Int_t iHt=0; iHt<5; iHt++){
	  //cfd
	  if(allData[ik+1][iHt]>0)
	    GetRawsData(shift+ik+1) -> 
	      Fill(allData[ik+1][iHt]-allData[1][0]);
	  //led
	  if(allData[ik+13][iHt] > 0 && allData[13][iHt]>0) 
	    GetRawsData(shift+ik+24+1)->
	      Fill(allData[ik+13][iHt]-allData[13][iHt]);
	  //led -cfd

	  if(allData[ik+13][iHt] > 0 && allData[ik+1][iHt] >0 )
	    GetRawsData(shift+ik+48+1)->
	      Fill(allData[ik+13][iHt]-allData[ik+1][iHt]);
	  //qtc
	  if(allData[2*ik+25][iHt] > 0 && allData[2*ik+26][iHt] > 0)
	    GetRawsData(shift+ik+72+1)->
	      Fill(allData[2*ik+25][iHt]-allData[2*ik+26][iHt]);
	  
	  
	  if(type == 8 && allData[ik+1][iHt]>0 ) feffC[ik]++; 
	  if(type == 8 && allData[ik+13][iHt]>0 ) feffA[ik]++; 
	  // GetRawsData(206)->Fill(ik,allData[ik+13][iHt]-allData[13][0]);
	  if(type == 8 && (allData[2*ik+25][iHt]>0 && allData[2*ik+26][iHt]>0) ) 
	    GetRawsData(207)->Fill(ik,allData[2*ik+25][iHt]-allData[2*ik+26][iHt]);
	}
	effic = Float_t(feffC[ik])/Float_t(fnEvent);
	GetRawsData(205)->Fill(ik,effic );
	effic = Float_t(feffA[ik])/Float_t(fnEvent);
	GetRawsData(206)->Fill(ik,effic );
	//      printf("CFD  efficiency for PMT %i = %f \n", ik, effic); 
      }
      for (Int_t ik = 12; ik<24; ik++) {
	for (Int_t iHt=0; iHt<5; iHt++) {
	  if(allData[ik+45][iHt]>0) 
	    //cfd
	    GetRawsData(shift+ik+1)->
	      Fill(allData[ik+45][iHt]-allData[57][0]);
	    //qtc	  
	  if(allData[2*ik+57][iHt]>0 && allData[2*ik+58][iHt]>0)
	    GetRawsData(shift+ik+72+1)->
	      Fill(allData[2*ik+57][iHt]-allData[2*ik+58][iHt]);
	  //led
	  if(allData[ik+57][iHt] > 0 )
	    GetRawsData(shift+ik+24+1)->
	      Fill(allData[ik+57][iHt]-allData[69][iHt]);
	    //led-cfd
	  if(allData[ik+57][iHt] > 0 &&allData[ik+45][iHt]>0)
	    GetRawsData(shift+ik+48+1)->
	      Fill(allData[ik+57][iHt]-allData[ik+45][iHt]);
	  
	  if(type == 8 && allData[ik+45][iHt]>0 )  feffC[ik]++; 

	  if(type == 8 && allData[ik+57][iHt]>0 )  feffA[ik]++; 
	  // GetRawsData(206)->Fill(ik,allData[ik+57][iHt]-allData[69][0]);
	  if(type == 8 && (allData[2*ik+57][iHt]>0 && allData[2*ik+58][iHt]>0) )
	    GetRawsData(207)->Fill(ik,allData[2*ik+57][iHt]-allData[2*ik+58][iHt]);	  	  
	}
	effic = Float_t(feffC[ik])/Float_t(fnEvent);
	GetRawsData(205)->Fill(ik,effic );
	effic = Float_t(feffA[ik])/Float_t(fnEvent);
	GetRawsData(206)->Fill(ik,effic );
	//   printf("CFD  efficiency for PMT %i = %f \n", ik, effic);       }
      }
	
      Int_t trChannel[6] = {49,50,51,52,55,56};  
      if(type == 7)
	{
	  for (Int_t iHt=0; iHt<6; iHt++) {
	    for (Int_t itr=0; itr<6; itr++) {
	      if(allData[trChannel[itr]][iHt]>0) fNumTriggers[itr]++;
	    }
	  }
	}
      if(type == 8)
	{
	  for (Int_t iHt=0; iHt<5; iHt++) {
	    for (Int_t itr=0; itr<6; itr++) {
	      if(allData[trChannel[itr]][iHt]>0)
		{
		  
		  GetRawsData(198+itr)->Fill(allData[trChannel[itr]][iHt]-allData[1][0]);
		  
		  fNumTriggersCal[itr]++;
		}
	    }
	    if(allData[53][iHt]>0 && allData[54][iHt]>0) 
	      GetRawsData(204)->Fill(allData[53][iHt]-allData[54][iHt]);
	  }
	} 

      
      
      delete start;
      }
    }
  


//____________________________________________________________________________
void AliT0QADataMakerRec::MakeDigits( TTree *digitsTree)
{
  //fills QA histos for Digits
  
  TArrayI *digCFD = new TArrayI(24);
  TArrayI *digLED = new TArrayI(24);
  TArrayI *digQT0 = new TArrayI(24);
  TArrayI *digQT1 = new TArrayI(24);
  Int_t refpoint=0;
  
  TBranch *brDigits=digitsTree->GetBranch("T0");
  AliT0digit *fDigits = new AliT0digit() ;
  if (brDigits) {
    brDigits->SetAddress(&fDigits);
  }else{
    AliError(Form("EXEC Branch T0 digits not found"));
    return;
  }
  digitsTree->GetEvent(0);
  digitsTree->GetEntry(0);
  brDigits->GetEntry(0);
  fDigits->GetTimeCFD(*digCFD);
  fDigits->GetTimeLED(*digLED);
  fDigits->GetQT0(*digQT0);
  fDigits->GetQT1(*digQT1);
  refpoint = fDigits->RefPoint();
  for (Int_t i=0; i<24; i++)
    {
    if (digCFD->At(i)>0) {
      Int_t cfd=digCFD->At(i)- refpoint;
      GetDigitsData(0) ->Fill(i,cfd);
      GetDigitsData(1) -> Fill(i,(digLED->At(i) - digCFD->At(i)));
      GetDigitsData(2) -> Fill(i, (digQT1->At(i) - digQT0->At(i)));
    }
    }  
  
  delete digCFD;
  delete digLED;
  delete digQT0;
  delete digQT1;
  
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
    if(i<12)
      GetRecPointsData(0) -> Fill(i, frecpoints -> GetTime(i) - frecpoints -> GetTime(0)); 
    if(i>11)
      GetRecPointsData(0) -> Fill(i,  frecpoints -> GetTime(i) - frecpoints -> GetTime(12)); 
    GetRecPointsData(1) -> Fill( i, frecpoints -> GetAmp(i) - frecpoints->AmpLED(i));
  }
  Double_t mmm=frecpoints->GetOnlineMean()- frecpoints->GetMeanTime();
   GetRecPointsData(2) ->Fill(mmm);
   // printf(" AliT0QADataMakerRec: diff mean  %f \n",mmm ); 
   // GetRecPointsData(73) ->Fill(frecpoints->GetMeanTime());
  


}

//____________________________________________________________________________
void AliT0QADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //fills QA histos for ESD

  GetESDsData(0) -> Fill(esd->GetT0());
  GetESDsData(1)-> Fill(esd->GetT0zVertex());

}
