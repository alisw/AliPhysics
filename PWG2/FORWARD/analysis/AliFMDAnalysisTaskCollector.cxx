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

/* $Id:$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Here some comments on what this task does                                //
 
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TInterpreter.h>
#include <TList.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>

#include "AliFMDAnalysisTaskCollector.h"
#include "AliAnalysisManager.h"
#include "AliESDFMD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliMultiplicity.h"
#include "AliESDVertex.h"
#include "AliFMDAnaParameters.h"
//#include "AliFMDGeometry.h"


ClassImp(AliFMDAnalysisTaskCollector)

//____________________________________________________________________
Double_t  AliFMDAnalysisTaskCollector::TripleLandau(Double_t *x, Double_t *par) {
  
  Double_t energy        = x[0];
  Double_t constant = par[0];
  Double_t mpv      = par[1];
  Double_t sigma    = par[2];
  Double_t alpha    = par[3];
  Double_t beta     = par[4];
 
  Double_t f = constant*(TMath::Landau(energy,mpv,sigma,kTRUE)+
			 alpha*TMath::Landau(energy,2*mpv+2*sigma*TMath::Log(2),2*sigma,kTRUE)+
			 beta*TMath::Landau(energy,3*mpv+3*sigma*TMath::Log(3),3*sigma,kTRUE) );
  
  return f;
}
//____________________________________________________________________

AliFMDAnalysisTaskCollector::AliFMDAnalysisTaskCollector()
: fDebug(0),
  fOutputList(0),
  fArray(0),
  fZvtxDist(0),
  fEvents(0),
  fEmptyEventsAside(0),
  fEmptyEventsCside(0)
{
  // Default constructor
  
 
}
//____________________________________________________________________
AliFMDAnalysisTaskCollector::AliFMDAnalysisTaskCollector(const char* name):
    AliAnalysisTaskSE(name),
    fDebug(0),
    fOutputList(0),
    fArray(0),
    fZvtxDist(0),
    fEvents(0),
    fEmptyEventsAside(0),
    fEmptyEventsCside(0)
{
  // Default constructor
  
  DefineOutput(1, TList::Class());
}
//____________________________________________________________________
void AliFMDAnalysisTaskCollector::UserCreateOutputObjects()
{
  // Create the output container
  printf("AnalysisTaskFMD::CreateOutPutData() \n");
  
  fOutputList = new TList();//(TList*)GetOutputData(0);
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  fArray     = new TObjArray();
  fArray->SetName("FMD");
  fArray->SetOwner();
  TH1F* hEdist = 0;
  TH1F* hEmptyEdist = 0;
  TH1F* hRingEdist = 0;
  for(Int_t nEta = 0; nEta < pars->GetNetaBins(); nEta++) {
    TObjArray* etaArray = new TObjArray();
    fArray->AddAtAndExpand(etaArray,nEta);
    for(Int_t det =1; det<=3;det++)
      {
	TObjArray* detArray = new TObjArray();
	detArray->SetName(Form("FMD%d",det));
	etaArray->AddAtAndExpand(detArray,det);
	Int_t nRings = (det==1 ? 1 : 2);
	for(Int_t ring = 0;ring<nRings;ring++)
	  {
	    Char_t ringChar = (ring == 0 ? 'I' : 'O');
	    hEdist = new TH1F(Form("FMD%d%c_etabin%d",det,ringChar,nEta),Form("FMD%d%c_etabin%d",det,ringChar,nEta),200,0,6);
	    hEdist->SetXTitle("#Delta E / E_{MIP}");
	    fOutputList->Add(hEdist);
	    detArray->AddAtAndExpand(hEdist,ring);
	    
	  } 
      }
    
  }
  
  for(Int_t det =1; det<=3;det++)
    {
      Int_t nRings = (det==1 ? 1 : 2);
      for(Int_t ring = 0;ring<nRings;ring++)
	{
	  Char_t ringChar = (ring == 0 ? 'I' : 'O');
	  hRingEdist = new TH1F(Form("ringFMD%d%c",det,ringChar),Form("ringFMD%d%c",det,ringChar),200,0,6);
	  hRingEdist->SetXTitle("#Delta E / E_{MIP}");
	  fOutputList->Add(hRingEdist);
	  hEmptyEdist = new TH1F(Form("emptyFMD%d%c",det,ringChar),Form("emptyFMD%d%c",det,ringChar),200,0,6);
	  hEmptyEdist->SetXTitle("#Delta E / E_{MIP}");
	  fOutputList->Add(hEmptyEdist);
	  
	}
    }
  fZvtxDist  = new TH1F("ZvtxDist","Vertex distribution",100,-30,30);
  fZvtxDist->SetXTitle("z vertex");
  fOutputList->Add(fZvtxDist);
}

//____________________________________________________________________
void AliFMDAnalysisTaskCollector::UserExec(Option_t */*option*/)
{
  
  AliESDEvent* esd = (AliESDEvent*)InputEvent();
  AliESD* old = esd->GetAliESDOld();
  if (old) {
    esd->CopyFromOldESD();
  }
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  TString triggers = esd->GetFiredTriggerClasses();
  //if(!triggers.IsNull()) return;
  //Bool_t trigger = pars->IsEventTriggered(esd);
  
  Bool_t physics = pars->IsEventTriggered(esd);
  //Bool_t empty   = pars->IsEventTriggered(esd,AliFMDAnaParameters::kEMPTY);
  Bool_t emptyAside = triggers.Contains("CINT1A-ABCE-NOPF-ALL");
  Bool_t emptyCside = triggers.Contains("CINT1C-ABCE-NOPF-ALL");
  
  //std::cout<<physics<<"   "<<empty<<std::endl;
  //if(!trigger)
  //  physics = kFALSE;
  Double_t vertex[3];
  
  Bool_t vtxStatus =   pars->GetVertex(esd,vertex);
  if(!vtxStatus)
    physics = kFALSE;
  
  if(physics) {
    fZvtxDist->Fill(vertex[2]);
    if(TMath::Abs(vertex[2]) > pars->GetVtxCutZ())
      physics = kFALSE;
  }
  
  AliESDFMD* fmd = esd->GetFMDData();
  if (!fmd) return;
  if(physics)
    fEvents++;
  else if(emptyAside) 
    fEmptyEventsAside++;
  else if(emptyCside) 
    fEmptyEventsCside++;
  
  if(!physics && !emptyAside && !emptyCside)
    return;
  TH1F* Edist = 0;
  TH1F* emptyEdist = 0;
  TH1F* ringEdist = 0;
  for(UShort_t det=1;det<=3;det++) {
    
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      emptyEdist = (TH1F*)fOutputList->FindObject(Form("emptyFMD%d%c",det,ring));
      ringEdist = (TH1F*)fOutputList->FindObject(Form("ringFMD%d%c",det,ring));
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      
      for(UShort_t sec =0; sec < nsec;  sec++)  {
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  
	  
	  Float_t mult = fmd->Multiplicity(det,ring,sec,strip);
	  if(mult == AliESDFMD::kInvalidMult || mult == 0) continue;
	  
	  Float_t eta = pars->GetEtaFromStrip(det,ring,sec,strip,vertex[2]);
	  
	  Int_t nEta  = pars->GetEtaBin(eta);
	  //	  std::cout<<det<<"  "<<ring<<"   "<<sec<<"    "<<strip<<"   "<<vertex[2]<<"   "<<eta<<"   "<<nEta<<std::endl;
	  if(physics) {
	    Edist = (TH1F*)fOutputList->FindObject(Form("FMD%d%c_etabin%d",det,ring,nEta));	  
	    Edist->Fill(mult);
	    ringEdist->Fill(mult);
	  }
	  if(emptyAside && det == 3 /*&& ring == 'O'*/) {
	    emptyEdist->Fill(mult);
	  }
	  //if(emptyCside && det == 3 && ring == 'I') {
	  //  emptyEdist->Fill(mult);
	  //}
	  
	  

	  if((emptyAside || emptyCside) && !(det == 3 /*&& ring == 'O'*/)/* && !(det == 2 && ring == 'O')*/) {
	    emptyEdist->Fill(mult);
	  }
	  //if(emptyCside && det == 2 && ring == 'O') {
	  //  emptyEdist->Fill(mult);
	  // }
	 
		  
	  //else {
	  //  AliWarning("Something is wrong - wrong trigger");
	  //  continue;
	  //}
	 
	  
	}
      }
    }
    
    PostData(1, fOutputList); 
    
  }
}
//____________________________________________________________________
  
void AliFMDAnalysisTaskCollector::Terminate(Option_t */*option*/) {
  std::cout<<"Analysed "<<fEvents<<" events and "<<fEmptyEventsAside<<" empty from A side and "<<fEmptyEventsCside<<"  on the C side"<<std::endl;
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
 
    
    for(Int_t det = 1; det<=3; det++) {
      Int_t nRings  =  (det == 1 ? 1 : 2);
      for(Int_t ring = 0;ring<nRings; ring++) {
	
	Char_t ringChar = (ring == 0 ? 'I' : 'O');
	TH1F* hRingEdist = (TH1F*)fOutputList->FindObject(Form("ringFMD%d%c",det,ringChar));
	TH1F* hEmptyEdist = (TH1F*)fOutputList->FindObject(Form("emptyFMD%d%c",det,ringChar));
	if(fEmptyEventsAside && det == 3 /*&& ringChar == 'O'*/)
	  hEmptyEdist->Scale(1./(Float_t)fEmptyEventsAside);
	
	//if(fEmptyEventsAside && det == 3 && ringChar == 'I')
	//  hEmptyEdist->Scale(1./(Float_t)fEmptyEventsCside);
	//if(fEmptyEventsAside && det == 2 && ringChar == 'O') //A side
	//  hEmptyEdist->Scale(1./(Float_t)fEmptyEventsAside);
	
	if((fEmptyEventsAside != 0 || fEmptyEventsCside !=0) && !(det == 3/* && ringChar == 'O'*/)/* && !(det == 2 && ringChar == 'O')*/)
	  hEmptyEdist->Scale(1./((Float_t)fEmptyEventsAside + ((Float_t)fEmptyEventsCside)));
	
	if(fEvents)
	  hRingEdist->Scale(1./(Float_t)fEvents);
	for(Int_t nEta = 0; nEta < pars->GetNetaBins(); nEta++) {
	  TH1F* hEdist = (TH1F*)fOutputList->FindObject(Form("FMD%d%c_etabin%d",det,ringChar,nEta));
	  if(fEvents)
	    hEdist->Scale(1./(Float_t)fEvents);
	  
	
      }
    }
  }
  
}
//____________________________________________________________________
void AliFMDAnalysisTaskCollector::ReadFromFile(const Char_t* filename, Bool_t store, Int_t speciesOption)  {

  //speciesOption:
  //0: p+p Landau fit
  //1: Pb+Pb triple landau convolution fit
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  pars->Init(kTRUE,AliFMDAnaParameters::kBackgroundCorrection);
  
  TFile fin(filename,"UPDATE");
  
  TList* list = (TList*)fin.Get("energyDist");
  
  AliFMDAnaCalibEnergyDistribution* EnergyDist = new AliFMDAnaCalibEnergyDistribution();
  
  EnergyDist->SetNetaBins(pars->GetNetaBins());
  EnergyDist->SetEtaLimits(pars->GetEtaMin(),pars->GetEtaMax());
    
  TF1* fitFunc = 0;
  
  for(Int_t det = 1; det<=3; det++) {
    Int_t nRings  =  (det == 1 ? 1 : 2);
    for(Int_t ring = 0;ring<nRings; ring++) {
      Char_t ringChar = (ring == 0 ? 'I' : 'O');
      TH1F* hEmptyEdist = (TH1F*)list->FindObject(Form("emptyFMD%d%c",det,ringChar));
      TH1F* hRingEdist = (TH1F*)list->FindObject(Form("ringFMD%d%c",det,ringChar));
      fitFunc = FitEnergyDistribution(hEmptyEdist, speciesOption) ;
      if(fitFunc)
	fitFunc->Write(Form("emptyFMD%d%c_fitfunc",det,ringChar),TObject::kWriteDelete);
      fitFunc = FitEnergyDistribution(hRingEdist, speciesOption) ;
      if(fitFunc)
	fitFunc->Write(Form("FMD%d%c_fitfunc",det,ringChar),TObject::kWriteDelete);
      
      
      EnergyDist->SetEmptyEnergyDistribution(det,ringChar,hEmptyEdist);
      EnergyDist->SetRingEnergyDistribution(det,ringChar,hRingEdist);
      for(Int_t nEta = 0; nEta < pars->GetNetaBins(); nEta++) {
	TH1F* hEdist = (TH1F*)list->FindObject(Form("FMD%d%c_etabin%d",det,ringChar,nEta));
	
	fitFunc = FitEnergyDistribution(hEdist, speciesOption) ;
	if(fitFunc)
	  fitFunc->Write(Form("FMD%d%c_etabin%d_fitfunc",det,ringChar,nEta),TObject::kWriteDelete);
	EnergyDist->SetEnergyDistribution(det,ringChar,nEta,hEdist);
	
      }
      
    }
    
  }
  
  fin.Write();
  fin.Close();
  
  if(store) {
    TFile fcalib(pars->GetPath(pars->GetEdistID() ),"RECREATE");
    EnergyDist->Write(AliFMDAnaParameters::GetEdistID());
    fcalib.Close();
  }



}
//____________________________________________________________________
TF1* AliFMDAnalysisTaskCollector::FitEnergyDistribution(TH1F* hEnergy, Int_t speciesOption) {
  
  TF1* fitFunc = 0;
  if(hEnergy->GetEntries() != 0) {
	  
    hEnergy->GetXaxis()->SetRangeUser(0.2,hEnergy->GetXaxis()->GetXmax());
    
    if(speciesOption == 0)
      fitFunc =  new TF1("FMDfitFunc","landau",hEnergy->GetBinCenter(hEnergy->GetMaximumBin())-0.2,3);
    if(speciesOption == 1) {
      fitFunc = new TF1("FMDfitFunc",TripleLandau,hEnergy->GetBinCenter(hEnergy->GetMaximumBin())-0.2,5,5);
      fitFunc->SetParNames("constant","MPV","sigma","2-Mip weight","3-Mip weight");
      fitFunc->SetParameters(10,0.8,0.1,0.05,0.01);
      fitFunc->SetParLimits(1,0.6,1.2);
      fitFunc->SetParLimits(3,0,1);
      fitFunc->SetParLimits(4,0,1);
      
    }
    hEnergy->Fit(fitFunc,"","",hEnergy->GetBinCenter(hEnergy->GetMaximumBin())-0.2,3);
    
  }
  return fitFunc;
  
}
//____________________________________________________________________
//
// EOF
//
