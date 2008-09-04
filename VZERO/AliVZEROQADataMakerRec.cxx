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


/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
*/

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2I.h> 
#include <TParameter.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliVZEROQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliVZERORawStream.h"
#include "AliVZEROReconstructor.h"


ClassImp(AliVZEROQADataMakerRec)
           
//____________________________________________________________________________ 
  AliVZEROQADataMakerRec::AliVZEROQADataMakerRec() : 
    AliQADataMakerRec(AliQA::GetDetName(AliQA::kVZERO), "VZERO Quality Assurance Data Maker"),
    fCalibData(GetCalibData()),
    fEvent(0)
    
{
  // constructor
  
   for(Int_t i=0; i<64; i++){  
       fEven[i] = 0;   
       fOdd[i]  = 0;  }
  
   for(Int_t i=0; i<128; i++){  
       fADC_Mean[i] = 0.0;   }	
}

//____________________________________________________________________________ 
  AliVZEROQADataMakerRec::AliVZEROQADataMakerRec(const AliVZEROQADataMakerRec& qadm) :
  AliQADataMakerRec(),
    fCalibData(GetCalibData()),
    fEvent(0)
  
{
  //copy constructor 
   SetName((const char*)qadm.GetName()) ; 
   SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliVZEROQADataMakerRec& AliVZEROQADataMakerRec::operator = (const AliVZEROQADataMakerRec& qadm )
{
  // Equal operator
  this->~AliVZEROQADataMakerRec();
  new(this) AliVZEROQADataMakerRec(qadm);
  return *this;
}

//____________________________________________________________________________
AliVZEROCalibData* AliVZEROQADataMakerRec::GetCalibData() const

{
  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  entry = man->Get("VZERO/Calib/Data");

  // Retrieval of data in directory VZERO/Calib/Data:

  AliVZEROCalibData *calibdata = 0;

  if (entry) calibdata = (AliVZEROCalibData*) entry->GetObject();
  if (!calibdata)  AliFatal("No calibration data from calibration database !");

  return calibdata;
}
 
//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kVZERO, task, list) ;
}

//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::InitESDs()
{
  //Create histograms to control ESD
 
  TH1I * h1 = new TH1I("hVZERONbPMA", "Number of PMs fired in V0A", 80, 0, 80); 
  h1->Sumw2() ;
  Add2ESDsList(h1, 0)  ;  
                                                                                                        
  TH1I * h2 = new TH1I("hVZERONbPMC", "Number of PMs fired in V0C", 80, 0, 80); 
  h2->Sumw2() ;
  Add2ESDsList(h2, 1) ;
 
  TH1I * h3 = new TH1I("hVZEROMultA", "Multiplicity in V0A", 50, 0, 50) ; 
  h3->Sumw2() ;
  Add2ESDsList(h3, 2) ;
 
  TH1I * h4 = new TH1I("hVZEROMultC", "Multiplicity in V0C", 50, 0, 50) ; 
  h4->Sumw2() ;
  Add2ESDsList(h4, 3) ;

  TH2F * h5 = new TH2F("fVzeroMult", "Vzero multiplicity",
			 64, -0.5, 63.5,1000, -0.5, 99.5);
  h5->GetXaxis()->SetTitle("Vzero PMT");
  h5->GetYaxis()->SetTitle("Multiplicity");
  h5->Sumw2() ;
  Add2ESDsList(h5, 4) ;
  TH1F * h6 = new TH1F("fBBA","BB Vzero A", 32, -0.5,31.5);
  h6->Sumw2();
  Add2ESDsList(h6, 5) ;
  TH1F * h7 = new TH1F("fBGA","BG Vzero A", 32, -0.5,31.5);
  h7->Sumw2();
  Add2ESDsList(h7, 6) ;
  TH1F * h8 = new TH1F("fBBC","BB Vzero C", 32, -0.5,31.5);
  h8->Sumw2();
  Add2ESDsList(h8, 7) ;
  TH1F * h9 = new TH1F("fBGC","BG Vzero C", 32, -0.5,31.5);
  h9->Sumw2();
  Add2ESDsList(h9, 8) ;

  TH2F *h10 = new TH2F("fVzeroAdc", "Vzero Adc",
			 64, -0.5, 63.5,1024, -0.5, 1023.5);
  h10->GetXaxis()->SetTitle("Vzero PMT");
  h10->GetYaxis()->SetTitle("Adc counts");
  h10->Sumw2() ;
  Add2ESDsList(h10, 9);

  TH2F *h11 = new TH2F("fVzeroTime", "Vzero Time",
			 64, -0.5, 63.5,300, -0.5, 149.5);
  h11->GetXaxis()->SetTitle("Vzero PMT");
  h11->GetYaxis()->SetTitle("Time [100 ps]");
  h11->Sumw2() ;
  Add2ESDsList(h11,10);

}

//____________________________________________________________________________ 
 void AliVZEROQADataMakerRec::InitRaws()
 {
   // create Raws histograms in Raws subdir

  char ADCname[12]; 
  char texte[40]; 
  TH1I *hRawADC0[64]; 
  TH1I *hRawADC1[64];
  
  Bool_t expert   = kTRUE ; 
  Bool_t saveCorr = kTRUE ; 
    
  TH2I *h0 = new TH2I("hCellADCMap0","ADC vs Cell for EVEN Integrator", 70, 0, 70, 512, 0, 1024);
  h0->Sumw2(); 
  Add2RawsList(h0,0, !expert, !saveCorr) ;
  TH2I *h1 = new TH2I("hCellADCMap1","ADC vs Cell for ODD Integrator", 70, 0, 70, 512, 0, 1024);
  h1->Sumw2();
  Add2RawsList(h1,1, !expert, !saveCorr) ;

  TH2F *hMeanADC0 = new TH2F("hCellMeanADCMap0","Mean ADC vs cell for EVEN integrator",70,-0.5,69.5,512,-0.5,511.5);       
  Add2RawsList(hMeanADC0,2, !expert, !saveCorr); 
  TH2F *hMeanADC1 = new TH2F("hCellMeanADCMap1","Mean ADC vs cell for ODD integrator",70,-0.5,69.5,512,-0.5,511.5);       
  Add2RawsList(hMeanADC1,3, !expert, !saveCorr);  
    
  TH1I *hMulV0A = new TH1I("hMulV0A", "Multiplicity in V0A", 40, 0, 40) ;  
  Add2RawsList(hMulV0A,4, !expert, saveCorr);  
  TH1I *hMulV0C = new TH1I("hMulV0C", "Multiplicity in V0C", 40, 0, 40) ;  
  Add2RawsList(hMulV0C,5, !expert, saveCorr);  
                          
  for (Int_t i=0; i<64; i++)
    {
       sprintf(ADCname,"hRaw0ADC%d",i);
       sprintf(texte,"Raw ADC in cell %d for even integrator",i);
       hRawADC0[i]= new TH1I(ADCname,texte,1024,0,1024);       
       Add2RawsList(hRawADC0[i],i+6, expert, !saveCorr);
                    
       sprintf(ADCname,"hRaw1ADC%d",i);
       sprintf(texte,"Raw ADC in cell %d for odd integrator",i);
       hRawADC1[i]= new TH1I(ADCname,texte,1024,0,1024);       
       Add2RawsList(hRawADC1[i],i+6+64, expert, !saveCorr);                     
     }  
 }

//____________________________________________________________________________
void AliVZEROQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs

  AliESDVZERO *esdVZERO=esd->GetVZEROData();
   
  if (esdVZERO) { 
      GetESDsData(0)->Fill(esdVZERO->GetNbPMV0A());
      GetESDsData(1)->Fill(esdVZERO->GetNbPMV0C());  
      GetESDsData(2)->Fill(esdVZERO->GetMTotV0A());
      GetESDsData(3)->Fill(esdVZERO->GetMTotV0C());  
      for(Int_t i=0;i<64;i++) {
	 GetESDsData(4)->Fill((Float_t) i,(Float_t) esdVZERO->GetMultiplicity(i));
	 GetESDsData(9)->Fill((Float_t) i,(Float_t) esdVZERO->GetAdc(i));
	 GetESDsData(10)->Fill((Float_t) i,(Float_t) esdVZERO->GetTime(i));
      }
      for(Int_t i=0;i<32;i++) { 
	if(esdVZERO->BBTriggerV0A(i)) 
	   GetESDsData(5)->Fill((Float_t) i);
	if(esdVZERO->BGTriggerV0A(i)) 
	   GetESDsData(6)->Fill((Float_t) i);
	if(esdVZERO->BBTriggerV0C(i)) 
	   GetESDsData(7)->Fill((Float_t) i);
	if(esdVZERO->BGTriggerV0C(i)) 
	   GetESDsData(8)->Fill((Float_t) i);
      }
  }
  
}

//____________________________________________________________________________
 void AliVZEROQADataMakerRec::MakeRaws(AliRawReader* rawReader)
 {
  //Fill histograms with Raws, computes average ADC values dynamically (pedestal subtracted)
                  
  rawReader->Reset() ; 
  AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader); 
  rawStream->Next();
  
  Float_t ChargeEoI, Threshold;  // for pedestal subtraction 
                              
//  for(Int_t i=0; i<128; i++) { printf(" i = %d pedestal = %f sigma = %f \n\n", 
//                                        i,  fCalibData->GetPedestal(i),fCalibData->GetSigma(i) );} 
             
  GetRawsData(2)->Reset();      // to keep only the last value of the ADC average
  GetRawsData(3)->Reset();

  Int_t mulV0A = 0 ; 
  Int_t mulV0C = 0 ; 
   
  Int_t thresholV0A = 50 ; 
  Int_t thresholV0C = 50 ; 
   
  for(Int_t j=0; j<64; j++) {
      Int_t i =  rawStream->GetOfflineChannel(j);  // Convert Online to Offline channel number 
      if(!rawStream->GetIntegratorFlag(i,10)) 
           { 
	    // even integrator - fills index6 to 69 
	    GetRawsData(0)->Fill(i,rawStream->GetADC(i)) ; 
      GetRawsData(i+6)->Fill(rawStream->GetADC(i)) ; 
	    ChargeEoI  = (float)rawStream->GetADC(i) - fCalibData->GetPedestal(i); 
	    Threshold  = 3.0 * fCalibData->GetSigma(i);   
	    if(rawStream->GetADC(i)<1023 && ChargeEoI > Threshold) {
	       fADC_Mean[i] = ((fADC_Mean[i]*fEven[i]) + rawStream->GetADC(i)) / (fEven[i]+1); 
	       GetRawsData(2)->Fill((Float_t) i, fADC_Mean[i]);
	       fEven[i]+=1; }
	   }
      else if(rawStream->GetIntegratorFlag(i,10)) 
           { 
	    // odd integrator  - fills index 70 to 133	
      GetRawsData(1)->Fill(i,rawStream->GetADC(i)) ; 
      GetRawsData(i+6+64)->Fill(rawStream->GetADC(i)) ;
	    ChargeEoI  = (float)rawStream->GetADC(i) - fCalibData->GetPedestal(i+64) ; 
	    Threshold  = 3.0 * fCalibData->GetSigma(i+64);  
	    if(rawStream->GetADC(i)<1023 && ChargeEoI > Threshold) { 
	       fADC_Mean[i+64] = ((fADC_Mean[i+64]*fOdd[i]) + rawStream->GetADC(i)) / (fOdd[i]+1);
	       GetRawsData(3)->Fill((Float_t) i, fADC_Mean[i+64]);
	       fOdd[i]+=1; }
	   } 	
     if( i < 32 ){
       if (rawStream->GetADC(i) > thresholV0C) 
           mulV0C++ ;     
     }else{
       if (rawStream->GetADC(i) > thresholV0A) 
           mulV0A++ ; 
     }
  }   
   
  fEvent++; 
  GetRawsData(4)->Fill(mulV0A) ; 
  TParameter<double> * p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(4)->GetName()))) ; 
  p->SetVal(mulV0A) ; 
  GetRawsData(5)->Fill(mulV0C) ; 
  p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(5)->GetName()))) ; 
  p->SetVal(mulV0C) ;                     
 }

//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::StartOfDetectorCycle()
{
  // Detector specific actions at start of cycle
  
}
