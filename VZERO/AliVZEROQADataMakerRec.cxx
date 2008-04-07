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

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliVZEROQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliVZERORawStream.h"
#include "AliVZEROReconstructor.h"


ClassImp(AliVZEROQADataMakerRec)
           
//____________________________________________________________________________ 
  AliVZEROQADataMakerRec::AliVZEROQADataMakerRec() : 
  AliQADataMakerRec(AliQA::GetDetName(AliQA::kVZERO), "VZERO Quality Assurance Data Maker")
{
  // constructor
}

//____________________________________________________________________________ 
AliVZEROQADataMakerRec::AliVZEROQADataMakerRec(const AliVZEROQADataMakerRec& qadm) :
  AliQADataMakerRec()
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
void AliVZEROQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kVZERO, task, list) ;  
}

//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::InitESDs()
{
  //Create histograms to controll ESD
 
  TH1I * h1 = new TH1I("hVZERONbPMA", "Number of PMs fired in V0A", 80, 0, 79); 
  h1->Sumw2() ;
  Add2ESDsList(h1, 0)  ;  
                                                                                                        
  TH1I * h2 = new TH1I("hVZERONbPMC", "Number of PMs fired in V0C", 80, 0, 79); 
  h2->Sumw2() ;
  Add2ESDsList(h2, 1) ;
 
  TH1I * h3 = new TH1I("hVZEROMultA", "Multiplicity in V0A", 50, 0, 49) ; 
  h3->Sumw2() ;
  Add2ESDsList(h3, 2) ;
 
  TH1I * h4 = new TH1I("hVZEROMultC", "Multiplicity in V0C", 50, 0, 49) ; 
  h4->Sumw2() ;
  Add2ESDsList(h4, 3) ;
	
}


//____________________________________________________________________________ 
 void AliVZEROQADataMakerRec::InitRaws()
 {
   // create Raws histograms in Raws subdir

  char ADCname[12]; 
  char texte[40]; 
  TH1I *fhRawADC0[64]; 
  TH1I *fhRawADC1[64]; 
  
  TH2I * h0 = new TH2I("hCellADCMap0","ADC vs Cell for EVEN Integrator", 70, 0, 69, 512, 0, 1023);
  h0->Sumw2(); 
  Add2RawsList(h0,0) ;
  TH2I * h1 = new TH2I("hCellADCMap1","ADC vs Cell for ODD Integrator", 70, 0, 69, 512, 0, 1023);
  h1->Sumw2();
  Add2RawsList(h1,1) ;
                           
  for (Int_t i=0; i<64; i++)
    {
       sprintf(ADCname,"hRaw0ADC%d",i);
       sprintf(texte,"Raw ADC in cell %d for even integrator",i);
       fhRawADC0[i]= new TH1I(ADCname,texte,1024,0,1023);       
       Add2RawsList(fhRawADC0[i],i+2);
       
       sprintf(ADCname,"hRaw1ADC%d",i);
       sprintf(texte,"Raw ADC in cell %d for odd integrator",i);
       fhRawADC1[i]= new TH1I(ADCname,texte,1024,0,1023);       
       Add2RawsList(fhRawADC1[i],i+2+64);       
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
  }

}

//____________________________________________________________________________
 void AliVZEROQADataMakerRec::MakeRaws(AliRawReader* rawReader)
 {
   //Fill histograms with Raws
       
  AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader); 
  rawStream->Next();
    
  for(Int_t i=0; i<64; i++) 
     {
       if(!rawStream->GetIntegratorFlag(i,10)) 
         { 
	   // even integrator - fills index 2 to 65   
	   GetRawsData(0)->Fill(i,rawStream->GetADC(i)) ; 
           GetRawsData(i+2)->Fill(rawStream->GetADC(i)) ; }
       else 
         { 
	   // odd integrator  - fills index 66 to 129	
           GetRawsData(1)->Fill(i,rawStream->GetADC(i)) ; 
           GetRawsData(i+2+64)->Fill(rawStream->GetADC(i)) ; }  	
      }          
 }

//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::StartOfDetectorCycle()
{
  // Detector specific actions at start of cycle
  
}
