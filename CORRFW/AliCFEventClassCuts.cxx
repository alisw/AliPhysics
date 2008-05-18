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
// To select on the Event Class: look at the Trigger mask and the ZDC info.
// Only pp-running trigger types implemented so far
// handles all masks for the trigger description
// and some general combinations like MB1,MB2,MB3,MB4 and MB5.
// The argument of IsSelected member function (passed object) is cast into 
// an AliVEvent, but cuts have a true meaning only for AliESD(AOD)Event 
// type objects.
// The class derives from AliCFCutBase
// Author:S.Arcelli Silvia.Arcelli@cern.ch
//
//
#include "TH1F.h"
#include "TList.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliCFEventClassCuts.h"
ClassImp(AliCFEventClassCuts) 
//____________________________________________________________________
AliCFEventClassCuts::AliCFEventClassCuts() : 
  AliCFCutBase(),
  fTriggerType(0),
  fTriggerAND(kFALSE),
  fZDCN1EnergyMin(-1.e99),  
  fZDCP1EnergyMin(-1.e99),  
  fZDCN2EnergyMin(-1.e99),  
  fZDCP2EnergyMin(-1.e99),  
  fZDCEM1EnergyMin(-1.e99),  
  fZDCEM2EnergyMin(-1.e99),  
  fZDCN1EnergyMax(1.e99),  
  fZDCP1EnergyMax(1.e99),  
  fZDCN2EnergyMax(1.e99),  
  fZDCP2EnergyMax(1.e99),  
  fZDCEM1EnergyMax(1.e99),  
  fZDCEM2EnergyMax(1.e99),
  fBitMap(0x0)
{
  //
  //ctor
  //  

  fBitMap=new TBits(0);
  Initialise();
}

//____________________________________________________________________
AliCFEventClassCuts::AliCFEventClassCuts(Char_t* name, Char_t* title) : 
  AliCFCutBase(name,title),
  fTriggerType(0),
  fTriggerAND(kFALSE),
  fZDCN1EnergyMin(-1.e99),  
  fZDCP1EnergyMin(-1.e99),  
  fZDCN2EnergyMin(-1.e99),  
  fZDCP2EnergyMin(-1.e99),  
  fZDCEM1EnergyMin(-1.e99),  
  fZDCEM2EnergyMin(-1.e99),  
  fZDCN1EnergyMax(1.e99),  
  fZDCP1EnergyMax(1.e99),  
  fZDCN2EnergyMax(1.e99),  
  fZDCP2EnergyMax(1.e99),  
  fZDCEM1EnergyMax(1.e99),  
  fZDCEM2EnergyMax(1.e99), 
  fBitMap(0x0)
{
  //
  //ctor
  //
  fBitMap=new TBits(0);
  Initialise();
 }

//_____________________________________________________________________________
AliCFEventClassCuts::AliCFEventClassCuts(const AliCFEventClassCuts& c) : 
  AliCFCutBase(c),
  fTriggerType(c.fTriggerType),
  fTriggerAND(c.fTriggerAND),
  fZDCN1EnergyMin(c.fZDCN1EnergyMin),  
  fZDCP1EnergyMin(c.fZDCP1EnergyMin),  
  fZDCN2EnergyMin(c.fZDCN2EnergyMin),  
  fZDCP2EnergyMin(c.fZDCP2EnergyMin),  
  fZDCEM1EnergyMin(c.fZDCEM1EnergyMin),  
  fZDCEM2EnergyMin(c.fZDCEM2EnergyMin),  
  fZDCN1EnergyMax(c.fZDCN1EnergyMax),  
  fZDCP1EnergyMax(c.fZDCP1EnergyMax),  
  fZDCN2EnergyMax(c.fZDCN2EnergyMax),  
  fZDCP2EnergyMax(c.fZDCP2EnergyMax),  
  fZDCEM1EnergyMax(c.fZDCEM1EnergyMax), 
  fZDCEM2EnergyMax(c.fZDCEM2EnergyMax),
  fBitMap(c.fBitMap)

{
  //
  //copy constructor
  //
}

//_____________________________________________________________________________
AliCFEventClassCuts& AliCFEventClassCuts::operator=(const AliCFEventClassCuts& c){
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fTriggerType = c.fTriggerType ;
    fTriggerAND = c.fTriggerAND ;
    fZDCN1EnergyMin = c.fZDCN1EnergyMin;  
    fZDCP1EnergyMin = c.fZDCP1EnergyMin;  
    fZDCN2EnergyMin = c.fZDCN2EnergyMin;  
    fZDCP2EnergyMin = c.fZDCP2EnergyMin;  
    fZDCEM1EnergyMin = c.fZDCEM1EnergyMin;  
    fZDCEM2EnergyMin = c.fZDCEM2EnergyMin;  
    fZDCN1EnergyMax = c.fZDCN1EnergyMax;  
    fZDCP1EnergyMax = c.fZDCP1EnergyMax;  
    fZDCN2EnergyMax = c.fZDCN2EnergyMax;  
    fZDCP2EnergyMax = c.fZDCP2EnergyMax;  
    fZDCEM1EnergyMax = c.fZDCEM1EnergyMax;  
    fZDCEM2EnergyMax = c.fZDCEM2EnergyMax;  
    fBitMap          = c.fBitMap;
  }


  for (Int_t i=0; i<c.kNCuts; i++){
    for (Int_t j=0; j<c.kNStepQA; j++){
      if(c.fhQA[i][j]) fhQA[i][j] = (TH1F*)c.fhQA[i][j]->Clone();
    }
  }

  return *this ;
}

//_____________________________________________________________________________
AliCFEventClassCuts::~AliCFEventClassCuts()
{
  //
  // destructor
  //
  for (Int_t i=0; i<kNCuts; i++){
    for (Int_t j=0; j<kNStepQA; j++){
      if(fhQA[i][j]) delete fhQA[i][j];
    }
  }

  if(fBitMap)delete fBitMap;

}

//_____________________________________________________________________________
void AliCFEventClassCuts::Initialise()
{
  //
  //initialization
  //
  // sets pointers to histos to zero
  for(Int_t i=0; i<kNCuts; i++){
    for(Int_t j =0; j<kNStepQA; j++){
      fhQA[i][j]=0x0;
    }
  }
}

//____________________________________________________________________
Bool_t AliCFEventClassCuts::IsSelected(TObject* obj) {
  //
  //Check if the requested cuts are passed
  //

  SelectionBitMap(obj);

  if (fIsQAOn) FillHistograms(obj,0);
  Bool_t isSelected = kTRUE;

  for (UInt_t icut=0; icut<fBitMap->GetNbits();icut++)
	if(!fBitMap->TestBitNumber(icut)) isSelected = kFALSE;

  if (!isSelected) return kFALSE ;
  if (fIsQAOn) FillHistograms(obj,1);
  return kTRUE;
}

//____________________________________________________________________
void AliCFEventClassCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins)
{
  //
  //setting x-axis bin limits of QA histogram fhQA[index] 
  // 
  for(Int_t i=0;i<kNStepQA;i++){
    if(!fhQA[index][i]){AliWarning("non-existing histogram!");
    return;
    }
    fhQA[index][i]->GetXaxis()->Set(nbins,bins);
  }
}
//____________________________________________________________________
void AliCFEventClassCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax)
{
  //
  //setting x-axis bins and range of QA histogram fhQA[index] 
  // 
  for(Int_t i=0;i<kNStepQA;i++){
    if(!fhQA[index][i]){AliWarning("non-existing histogram!");
    return;
    }
    fhQA[index][i]->GetXaxis()->Set(nbins,xmin,xmax);
  }
}
//____________________________________________________________________
void AliCFEventClassCuts::SelectionBitMap(TObject* obj) {
  //
  //cut on trigger type (just pp running trigger types implemented so far)
  //and on the energy observed in the ZDC. The argument is cast into 
  //an AliVEvent, but has true meaning only for AliESDEvent type objects.
  //Check if the requested cuts are passed and return a bitmap
  //

  for(Int_t j=0;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE);

  AliVEvent* esd = dynamic_cast<AliVEvent *>(obj);
  if (!esd ) return;


  //now start checking the cuts
  //first assume the event will be accepted: 
  for(Int_t j=0;j<kNCuts;j++)fBitMap->SetBitNumber(j,kTRUE);


  //Check the trigger:

  //look at the Trigger mask in current event
  TBits *triggerBitMap=new TBits(0);
  TriggerBitMap(esd,triggerBitMap); 
  //now compare to what was requested as a Trigger:  
  if(fTriggerType.GetNbits()>0)fBitMap->SetBitNumber(0,kFALSE); //trigger required, initialize to false
  for(Int_t j=0;j<kNTriggers+kNTriggersMB;j++){
    if(fTriggerType.TestBitNumber(j)){
      if(!fTriggerAND){
	if(triggerBitMap->TestBitNumber(j) == fTriggerType.TestBitNumber(j)){
	  fBitMap->SetBitNumber(0,kTRUE); 

	  break;// @least one requested bit fired, ok
	}
      }else{
	if(!triggerBitMap->TestBitNumber(j)){
	  break;
	}	
      }
    }
  }
  
  delete triggerBitMap;
  //Then, cut on the energy observed in the ZDC 
  
  if( esd->GetZDCN1Energy()<fZDCN1EnergyMin || esd->GetZDCN1Energy()>fZDCN1EnergyMax)fBitMap->SetBitNumber(1,kFALSE); 
  if( esd->GetZDCP1Energy()<fZDCP1EnergyMin || esd->GetZDCP1Energy()>fZDCP1EnergyMax)fBitMap->SetBitNumber(2,kFALSE); 
  if( esd->GetZDCN2Energy()<fZDCN2EnergyMin || esd->GetZDCN2Energy()>fZDCN2EnergyMax)fBitMap->SetBitNumber(3,kFALSE);
  if( esd->GetZDCP2Energy()<fZDCP2EnergyMin || esd->GetZDCP2Energy()>fZDCP2EnergyMax)fBitMap->SetBitNumber(4,kFALSE); 
  if( esd->GetZDCEMEnergy(0)<fZDCEM1EnergyMin || esd->GetZDCEMEnergy(0)>fZDCEM1EnergyMax)fBitMap->SetBitNumber(5,kFALSE); 
  if( esd->GetZDCEMEnergy(1)<fZDCEM2EnergyMin || esd->GetZDCEMEnergy(1)>fZDCEM2EnergyMax)fBitMap->SetBitNumber(6,kFALSE); 
  return ;

}

//_____________________________________________________________________________
Bool_t AliCFEventClassCuts::IsTriggered(AliVEvent* ev, TriggerType trigger) {
  //
  //look at the Trigger mask in current event
  TBits *triggerBitMap=new TBits(0);
  TriggerBitMap(ev,triggerBitMap); 
  Bool_t isTriggered=kFALSE;  
  if(triggerBitMap->TestBitNumber(trigger))isTriggered=kTRUE;
  delete triggerBitMap;
  return isTriggered;

}

//_____________________________________________________________________________
void AliCFEventClassCuts::TriggerBitMap(AliVEvent* ev, TBits *bitmapT ) {
  //

  for(Int_t itrig=0;itrig<kNTriggers+kNTriggersMB;itrig++)bitmapT->SetBitNumber(itrig,kFALSE);
  if (!ev ) return;

  ULong64_t triggerMask = ev->GetTriggerMask();
  //run over the different triggers in the mask, and check which bits have fired    
  for(Int_t itrig=0;itrig<kNTriggers;itrig++){
    bitmapT->SetBitNumber(itrig,kFALSE);
    if (triggerMask&(0x1 <<itrig)){
      bitmapT->SetBitNumber(itrig,kTRUE);
    }
  }

  //Trigger combinations, Minimum bias triggers

  //MB1 case: (GFO || V0OR) && !BG
  if((bitmapT->TestBitNumber(5) || (bitmapT->TestBitNumber(0) || bitmapT->TestBitNumber(1))) && !bitmapT->TestBitNumber(2)) bitmapT->SetBitNumber(17,kTRUE); 
 
  //MB2 case: (GFO && V0OR) && !BG
  if((bitmapT->TestBitNumber(5) && (bitmapT->TestBitNumber(0) || bitmapT->TestBitNumber(1))) && !bitmapT->TestBitNumber(2)) bitmapT->SetBitNumber(18,kTRUE); 

  //MB3 case : (GFO && V0AND) && !BG
  if((bitmapT->TestBitNumber(5) && (bitmapT->TestBitNumber(0) && bitmapT->TestBitNumber(1))) && !bitmapT->TestBitNumber(2)) bitmapT->SetBitNumber(19,kTRUE); 

  //MB4 case: (GFO || V0AND) && !BG
  if((bitmapT->TestBitNumber(5) || (bitmapT->TestBitNumber(0) && bitmapT->TestBitNumber(1))) &&  !bitmapT->TestBitNumber(2)) bitmapT->SetBitNumber(20,kTRUE); 

  //MB5 case:: GFO && !BG
  if(bitmapT->TestBitNumber(5) && !bitmapT->TestBitNumber(2)) bitmapT->SetBitNumber(21,kTRUE); 

  return;
} 

//_____________________________________________________________________________
 void AliCFEventClassCuts::DefineHistograms() {
  //
  // histograms for cut variables
  //
  Int_t color = 2;

  if(!fIsQAOn) {
    AliInfo(Form("No QA histos requested, Please first set the QA flag on!"));
    return;
  }  
  
  // book QA histograms
  Char_t str[256];
  for (Int_t i=0; i<kNStepQA; i++) {
    if (i==0) sprintf(str," ");
    else sprintf(str,"_cut");

    fhQA[kTrigger][i]	= new  TH1F(Form("%s_TriggerBits%s",GetName(),str),	                "",23,-0.5,22.5);
    fhQA[kZDCEnergyN1][i]	= new  TH1F(Form("%s_ZDC_Energy_N1%s",GetName(),str),		"",800,-500,7500);
    fhQA[kZDCEnergyP1][i]	= new  TH1F(Form("%s_ZDC_Energy_P1%s",GetName(),str),		"",800,-500,7500);
    fhQA[kZDCEnergyN2][i]	= new  TH1F(Form("%s_ZDC_Energy_N2%s",GetName(),str),		"",800,-500,7500);
    fhQA[kZDCEnergyP2][i]	= new  TH1F(Form("%s_ZDC_Energy_P2%s",GetName(),str),		"",800,-500,7500);
    fhQA[kZDCEnergyEM1][i]	= new  TH1F(Form("%s_ZDC_Energy_EM1%s",GetName(),str),		"",800,-500,7500);
    fhQA[kZDCEnergyEM2][i]	= new  TH1F(Form("%s_ZDC_Energy_EM2%s",GetName(),str),		"",800,-500,7500);

    fhQA[kTrigger][i]	        ->SetXTitle("Trigger Bits");
    fhQA[kZDCEnergyN1][i]	->SetXTitle("ZDC Energy N1 (GeV)");
    fhQA[kZDCEnergyP1][i]	->SetXTitle("ZDC Energy P1 (GeV)");
    fhQA[kZDCEnergyN2][i]	->SetXTitle("ZDC Energy N2 (GeV)");
    fhQA[kZDCEnergyP2][i]	->SetXTitle("ZDC Energy P2 (GeV)");
    fhQA[kZDCEnergyEM1][i]	->SetXTitle("ZDC Energy EM1 (GeV)");
    fhQA[kZDCEnergyEM2][i]	->SetXTitle("ZDC Energy EM2 (GeV)");

  }

  for(Int_t i=0; i<kNCuts; i++) fhQA[i][1]->SetLineColor(color);

}

//_____________________________________________________________________________
void AliCFEventClassCuts::FillHistograms(TObject* obj, Bool_t afterCuts)
{
  //
  // fill the QA histograms
  //
  if(!fIsQAOn) return;

  // cast TObject into VParticle
  AliVEvent* esd = dynamic_cast<AliVEvent *>(obj);
  if (!esd ) return  ;

  //look at the Trigger mask in current event
  TBits *triggerBitMap=new TBits(0);
  TriggerBitMap(esd, triggerBitMap); 
  
  //trigger Mask
  for(Int_t itrig=0;itrig<kNTriggers+kNTriggersMB;itrig++){
    if(triggerBitMap->TestBitNumber(itrig)){
      fhQA[kTrigger][afterCuts]->Fill(itrig);
    }
  }   

  delete triggerBitMap;

  //ZDC Quantities
  fhQA[kZDCEnergyN1][afterCuts] ->Fill(esd->GetZDCN1Energy());
  fhQA[kZDCEnergyP1][afterCuts] ->Fill(esd->GetZDCP1Energy());
  fhQA[kZDCEnergyN2][afterCuts] ->Fill(esd->GetZDCN2Energy());
  fhQA[kZDCEnergyP2][afterCuts] ->Fill(esd->GetZDCP2Energy());
  fhQA[kZDCEnergyEM1][afterCuts]->Fill(esd->GetZDCEMEnergy(0));
  fhQA[kZDCEnergyEM2][afterCuts]->Fill(esd->GetZDCEMEnergy(1));

}

//_____________________________________________________________________________
void AliCFEventClassCuts::AddQAHistograms(TList *qaList) {
  //
  // saves the histograms in a TList
  //

  DefineHistograms();

  for (Int_t j=0; j<kNStepQA; j++) {
    for(Int_t i=0; i<kNCuts; i++)
	qaList->Add(fhQA[i][j]);
  }
}
