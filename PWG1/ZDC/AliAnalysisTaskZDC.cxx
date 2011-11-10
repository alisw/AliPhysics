/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//							   //
//	Class to analyze ZDC data			   //
//							   //
/////////////////////////////////////////////////////////////

#include <TList.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TFile.h>
#include <TString.h>

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliMultiplicity.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskZDC.h"

ClassImp(AliAnalysisTaskZDC)


//________________________________________________________________________
AliAnalysisTaskZDC::AliAnalysisTaskZDC():
  AliAnalysisTaskSE(),
    fDebug(0),
    fIsMCInput(kFALSE),
    fOutput(0x0),
    fhTDCZNSum(0x0),
    fhTDCZNDiff(0x0),
    fhZNCSpectrum(0x0),
    fhZNASpectrum(0x0),
    fhZPCSpectrum(0x0),
    fhZPASpectrum(0x0),
    fhZEM1Spectrum(0x0),
    fhZEM2Spectrum(0x0),
    fhZNCpmc(0x0),	 
    fhZNApmc(0x0),	 
    fhZPCpmc(0x0),	 
    fhZPApmc(0x0),	 
    fhZNCCentroid(0x0), 
    fhZNACentroid(0x0), 
    fhZNCemd(0x0),	   
    fhZNAemd(0x0),
    fhPMCZNCemd(0x0), 
    fhPMCZNAemd(0x0),
    fDebunch(0x0)
{   
   // Default constructor
}   

//________________________________________________________________________
AliAnalysisTaskZDC::AliAnalysisTaskZDC(const char *name):
  AliAnalysisTaskSE(name),
    fDebug(0),
    fIsMCInput(kFALSE),
    fOutput(0x0),
    fhTDCZNSum(0x0),
    fhTDCZNDiff(0x0),
    fhZNCSpectrum(0x0),
    fhZNASpectrum(0x0),
    fhZPCSpectrum(0x0),
    fhZPASpectrum(0x0),
    fhZEM1Spectrum(0x0),
    fhZEM2Spectrum(0x0),
    fhZNCpmc(0x0),	 
    fhZNApmc(0x0),	 
    fhZPCpmc(0x0),	 
    fhZPApmc(0x0),	 
    fhZNCCentroid(0x0), 
    fhZNACentroid(0x0), 
    fhZNCemd(0x0),	   
    fhZNAemd(0x0),
    fhPMCZNCemd(0x0), 
    fhPMCZNAemd(0x0),
    fDebunch(0x0) 
{  
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class()); 

}

//________________________________________________________________________
AliAnalysisTaskZDC& AliAnalysisTaskZDC::operator=(const AliAnalysisTaskZDC& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c);
  }
  return *this;
}

//________________________________________________________________________
AliAnalysisTaskZDC::AliAnalysisTaskZDC(const AliAnalysisTaskZDC& ana):
  AliAnalysisTaskSE(ana),
  fDebug(ana.fDebug),	  
  fIsMCInput(ana.fIsMCInput),
  fOutput(ana.fOutput),
  fhTDCZNSum(ana.fhTDCZNSum),
  fhTDCZNDiff(ana.fhTDCZNDiff),
  fhZNCSpectrum(ana.fhZNCSpectrum),
  fhZNASpectrum(ana.fhZNASpectrum),
  fhZPCSpectrum(ana.fhZPCSpectrum),
  fhZPASpectrum(ana.fhZPASpectrum),
  fhZEM1Spectrum(ana.fhZEM1Spectrum),
  fhZEM2Spectrum(ana.fhZEM2Spectrum),
  fhZNCpmc(ana.fhZNCpmc),       
  fhZNApmc(ana.fhZNApmc),       
  fhZPCpmc(ana.fhZPCpmc),       
  fhZPApmc(ana.fhZPApmc),       
  fhZNCCentroid(ana.fhZNCCentroid), 
  fhZNACentroid(ana.fhZNACentroid), 
  fhZNCemd(ana.fhZNCemd),	 
  fhZNAemd(ana.fhZNAemd),
  fhPMCZNCemd(ana.fhPMCZNCemd), 
  fhPMCZNAemd(ana.fhPMCZNAemd),
  fDebunch(ana.fDebunch)
{
  //
  // Copy Constructor	
  //
}
 
//________________________________________________________________________
AliAnalysisTaskZDC::~AliAnalysisTaskZDC()
{
  // Destructor
  if(fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    delete fOutput; fOutput=0;
  } 
   
}  

//________________________________________________________________________
void AliAnalysisTaskZDC::UserCreateOutputObjects()
{
  // Create the output containers

  fOutput = new TList;
  fOutput->SetOwner();
  //fOutput->SetName("output");
  
  fhTDCZNSum = new TH1F("fhTDCZNSum","TDC_{ZNC}+TDC_{ZNA}",60,-30.,-30.);
  fhTDCZNSum->GetXaxis()->SetTitle("TDC_{ZNC}+TDC_{ZNA} (ns)");
  fOutput->Add(fhTDCZNSum);      
  
  fhTDCZNDiff = new TH1F("fhTDCZNDiff","TDC_{ZNC}-TDC_{ZNA}",60,-30.,30.);
  fhTDCZNDiff->GetXaxis()->SetTitle("TDC_{ZNC}-TDC_{ZNA} (ns)");
  fOutput->Add(fhTDCZNDiff);      
  
  fhZNCSpectrum = new TH1F("fhZNCSpectrum", "ZNC signal", 200,0., 140000.);
  fOutput->Add(fhZNCSpectrum);      
  fhZNASpectrum = new TH1F("fhZNASpectrum", "ZNA signal", 200,0., 140000.) ;
  fOutput->Add(fhZNASpectrum);      
  fhZPCSpectrum = new TH1F("fhZPCSpectrum", "ZPC signal", 200,0., 50000.) ;
  fOutput->Add(fhZPCSpectrum);      
  fhZPASpectrum = new TH1F("fhZPASpectrum", "ZPA signal", 200,0., 50000.) ;
  fOutput->Add(fhZPASpectrum);      
  fhZEM1Spectrum = new TH1F("fhZEM1Spectrum", "ZEM1 signal", 100,0., 2500.);
  fOutput->Add(fhZEM1Spectrum);      
  fhZEM2Spectrum = new TH1F("fhZEM2Spectrum", "ZEM2 signal", 100,0., 2500.);
  fOutput->Add(fhZEM2Spectrum);      
  
  fhZNCpmc = new TH1F("fhZNCpmc","ZNC PMC",200, 0., 160000.);
  fOutput->Add(fhZNCpmc);      
  fhZNApmc = new TH1F("fhZNApmc","ZNA PMC",200, 0., 160000.); 
  fOutput->Add(fhZNApmc);      
  fhZPCpmc = new TH1F("fhZPCpmc","ZPC PMC",200, 0., 40000.); 
  fOutput->Add(fhZPCpmc);      
  fhZPApmc = new TH1F("fhZPApmc","ZPA PMC",200, 0., 40000.); 
  fOutput->Add(fhZPApmc);      
  
  fhZNCCentroid = new TH2F("fhZNCCentroid","Centroid over ZNC",70,-3.5,3.5,70,-3.5,3.5); 
  fOutput->Add(fhZNCCentroid);      
  fhZNACentroid = new TH2F("fhZNACentroid","Centroid over ZNA",70,-3.5,3.5,70,-3.5,3.5); 
  fOutput->Add(fhZNACentroid);      
  
  fhZNCemd = new TH1F("fhZNCemd","ZNC signal lg",200,0.,6000.);	 
  fOutput->Add(fhZNCemd);      
  fhZNAemd = new TH1F("fhZNAemd","ZNA signal lg",200,0.,6000.);	 
  fOutput->Add(fhZNAemd);      
  fhPMCZNCemd = new TH1F("fhPMCZNCemd","ZNC PMC lg",200, 10., 6000.);   
  fOutput->Add(fhPMCZNCemd);      
  fhPMCZNAemd = new TH1F("fhPMCZNAemd","ZNA PMC lg",200, 10., 6000.);   
  fOutput->Add(fhPMCZNAemd);     
  
  fDebunch = new TH2F("fDebunch","ZN TDC sum vs. diff", 120,-30,30,120,-30,-30);
  fOutput->Add(fDebunch);     
    
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskZDC::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  if(fDebug>1) printf(" **** AliAnalysisTaskZDC::UserExec() \n");
  
  if (!InputEvent()) {
    Printf("ERROR: InputEvent not available");
    return;
  }

      
  AliESDEvent* esd = dynamic_cast<AliESDEvent*> (InputEvent());
  if(!esd) return;
  // Select PHYSICS events (type=7, for data)
  if(!fIsMCInput && esd->GetEventType()!=7) return; 
  
  // ********* MC INFO *********************************
  if(fIsMCInput){

    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }
  
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
   }

    AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
    if(!genHeader){
      printf("  Event generator header not available!!!\n");
      return;
    }

    /*if(genHeader->InheritsFrom(AliGenHijingEventHeader::Class())){
      Float_t bMC = ((AliGenHijingEventHeader*) genHeader)->ImpactParameter();
      Int_t specNeutronProj = ((AliGenHijingEventHeader*) genHeader)->ProjSpectatorsn();
      Int_t specProtonProj  = ((AliGenHijingEventHeader*) genHeader)->ProjSpectatorsp();
      Int_t specNeutronTarg = ((AliGenHijingEventHeader*) genHeader)->TargSpectatorsn();
      Int_t specProtonTarg  = ((AliGenHijingEventHeader*) genHeader)->TargSpectatorsp();
      Int_t npartTargMC = 208-(specNeutronTarg+specProtonTarg);
      Int_t npartProjMC = 208-(specNeutronProj+specProtonProj);
    }*/  

  }
  // ****************************************************
  
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    
  AliESDZDC *esdZDC = esd->GetESDZDC();
  
  if((((AliInputEventHandler*)(am->GetInputEventHandler()))->IsEventSelected())){
  
    fhZNCSpectrum->Fill(esdZDC->GetZDCN1Energy());	  
    fhZNASpectrum->Fill(esdZDC->GetZDCN2Energy());
    fhZPCSpectrum->Fill(esdZDC->GetZDCP1Energy());		  
    fhZPASpectrum->Fill(esdZDC->GetZDCP2Energy());	  
    fhZEM1Spectrum->Fill(esdZDC->GetZDCEMEnergy(0)/8.);
    fhZEM2Spectrum->Fill(esdZDC->GetZDCEMEnergy(1)/8.);
  
    const Double_t * towZNC = esdZDC->GetZN1TowerEnergy();
    const Double_t * towZPC = esdZDC->GetZP1TowerEnergy();
    const Double_t * towZNA = esdZDC->GetZN2TowerEnergy();
    const Double_t * towZPA = esdZDC->GetZP2TowerEnergy();
    //    
    fhZNCpmc->Fill(towZNC[0]);    
    fhZNApmc->Fill(towZNA[0]);    
    fhZPCpmc->Fill(towZPC[0]);    
    fhZPApmc->Fill(towZPA[0]);    
  
    Double_t xyZNC[2]={-99.,-99.}, xyZNA[2]={-99.,-99.};
    esdZDC->GetZNCentroidInPbPb(1380., xyZNC, xyZNA);
    //esdZDC->GetZNCentroidInpp(xyZNC, xyZNA);
    
    fhZNCCentroid->Fill(xyZNC[0], xyZNC[1]); 
    fhZNACentroid->Fill(xyZNA[0], xyZNA[1]); 
    
    const Double_t * towZNCLG = esdZDC->GetZN1TowerEnergyLR();
    const Double_t * towZNALG = esdZDC->GetZN2TowerEnergyLR();
    Double_t znclg=0., znalg=0.;
    for(Int_t iq=0; iq<5; iq++){
       znclg += towZNCLG[iq];
       znalg += towZNALG[iq];
    }
    fhZNCemd->Fill(znclg/2.);	 
    fhZNAemd->Fill(znalg/2.);	 
    fhPMCZNCemd->Fill(towZNCLG[0]);   
    fhPMCZNAemd->Fill(towZNALG[0]);   
  
  }
    
  Float_t tdcC=999., tdcA=999;
  Float_t tdcSum=999., tdcDiff=999;
  if(esdZDC->GetZDCTDCData(10,0)>1e-4){
    tdcC = esdZDC->GetZDCTDCCorrected(10,0)-esdZDC->GetZDCTDCCorrected(15,0);
    if(esdZDC->GetZDCTDCData(12,0)>1e-4){
      tdcA = esdZDC->GetZDCTDCCorrected(12,0)-esdZDC->GetZDCTDCCorrected(15,0);
      tdcSum = tdcC+tdcA;
      tdcDiff = tdcC-tdcA;
    }
  }
  //for(Int_t i=0; i<4; i++){
    if(tdcSum!=999.){
      fhTDCZNSum->Fill(tdcSum); 
      fDebunch->Fill(tdcDiff, tdcSum);  
    }
    if(tdcDiff!=999.)fhTDCZNDiff->Fill(tdcDiff); 
  //}
  
  PostData(1, fOutput);
   
}



//________________________________________________________________________
void AliAnalysisTaskZDC::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
}
