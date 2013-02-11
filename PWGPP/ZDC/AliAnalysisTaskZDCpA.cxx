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
#include "AliAnalysisTaskZDCpA.h"

ClassImp(AliAnalysisTaskZDCpA)


//________________________________________________________________________
AliAnalysisTaskZDCpA::AliAnalysisTaskZDCpA():
  AliAnalysisTaskSE(),
    fDebug(0),
    fIsMCInput(kFALSE),
    fOutput(0x0),
    fhTDCZNC(0x0),
    fhTDCZNA(0x0),
    fhTDCZNSum(0x0),
    fhTDCZNDiff(0x0),
    fhZNCSumQ(0x0),
    fhZNASumQ(0x0),
    fhZPCSumQ(0x0),
    fhZPASumQ(0x0),
    fhZEM1Spectrum(0x0),
    fhZEM2Spectrum(0x0),
    fhZNCpmc(0x0),	 
    fhZNApmc(0x0),	 
    fhZPCpmc(0x0),	 
    fhZPApmc(0x0),	 
    fhZNCCentroid(0x0), 
    fhZNACentroid(0x0), 
    fhPMCZNCemd(0x0), 
    fhPMCZNAemd(0x0),
    fDebunch(0x0),
    fhTDCZNAcorr(0x0),
    fhTDCZNCcorr(0x0)
{   
   // Default constructor
}   

//________________________________________________________________________
AliAnalysisTaskZDCpA::AliAnalysisTaskZDCpA(const char *name):
  AliAnalysisTaskSE(name),
    fDebug(0),
    fIsMCInput(kFALSE),
    fOutput(0x0),
    fhTDCZNC(0x0),
    fhTDCZNA(0x0),
    fhTDCZNSum(0x0),
    fhTDCZNDiff(0x0),
    fhZNCSumQ(0x0),
    fhZNASumQ(0x0),
    fhZPCSumQ(0x0),
    fhZPASumQ(0x0),
    fhZEM1Spectrum(0x0),
    fhZEM2Spectrum(0x0),
    fhZNCpmc(0x0),	 
    fhZNApmc(0x0),	 
    fhZPCpmc(0x0),	 
    fhZPApmc(0x0),	 
    fhZNCCentroid(0x0), 
    fhZNACentroid(0x0), 
    fhPMCZNCemd(0x0), 
    fhPMCZNAemd(0x0),
    fDebunch(0x0),
    fhTDCZNAcorr(0x0),
    fhTDCZNCcorr(0x0)

{  
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class()); 

}

//________________________________________________________________________
AliAnalysisTaskZDCpA& AliAnalysisTaskZDCpA::operator=(const AliAnalysisTaskZDCpA& c)
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
AliAnalysisTaskZDCpA::AliAnalysisTaskZDCpA(const AliAnalysisTaskZDCpA& ana):
  AliAnalysisTaskSE(ana),
  fDebug(ana.fDebug),	  
  fIsMCInput(ana.fIsMCInput),
  fOutput(ana.fOutput),
  fhTDCZNC(ana.fhTDCZNC),
  fhTDCZNA(ana.fhTDCZNA),
  fhTDCZNSum(ana.fhTDCZNSum),
  fhTDCZNDiff(ana.fhTDCZNDiff),
  fhZNCSumQ(ana.fhZNCSumQ),
  fhZNASumQ(ana.fhZNASumQ),
  fhZPCSumQ(ana.fhZPCSumQ),
  fhZPASumQ(ana.fhZPASumQ),
  fhZEM1Spectrum(ana.fhZEM1Spectrum),
  fhZEM2Spectrum(ana.fhZEM2Spectrum),
  fhZNCpmc(ana.fhZNCpmc),       
  fhZNApmc(ana.fhZNApmc),       
  fhZPCpmc(ana.fhZPCpmc),       
  fhZPApmc(ana.fhZPApmc),       
  fhZNCCentroid(ana.fhZNCCentroid), 
  fhZNACentroid(ana.fhZNACentroid), 
  fhPMCZNCemd(ana.fhPMCZNCemd), 
  fhPMCZNAemd(ana.fhPMCZNAemd),
  fDebunch(ana.fDebunch),
  fhTDCZNAcorr(ana.fhTDCZNAcorr),
  fhTDCZNCcorr(ana.fhTDCZNCcorr)

{
  //
  // Copy Constructor	
  //
}
 
//________________________________________________________________________
AliAnalysisTaskZDCpA::~AliAnalysisTaskZDCpA()
{
  // Destructor
  if(fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    delete fOutput; fOutput=0;
  } 
   
}  

//________________________________________________________________________
void AliAnalysisTaskZDCpA::UserCreateOutputObjects()
{
  // Create the output containers

  fOutput = new TList;
  fOutput->SetOwner();
  //fOutput->SetName("output");
  
  fhTDCZNC = new TH1F("fhTDCZNC","TDC_{ZNC}",120, -60., -20.);
  fhTDCZNC->GetXaxis()->SetTitle("TDC_{ZNC} (ns)");
  fOutput->Add(fhTDCZNC);      
  
  fhTDCZNA = new TH1F("fhTDCZNA","TDC_{ZNA}",120, -60., -20.);
  fhTDCZNA->GetXaxis()->SetTitle("TDC_{ZNA} (ns)");
  fOutput->Add(fhTDCZNA);     
  
  fhTDCZNSum = new TH1F("fhTDCZNSum","TDC_{ZNC}+TDC_{ZNA}",120,-30.,30.);
  fhTDCZNSum->GetXaxis()->SetTitle("TDC_{ZNC}+TDC_{ZNA} (ns)");
  fOutput->Add(fhTDCZNSum);      
  
  fhTDCZNDiff = new TH1F("fhTDCZNDiff","TDC_{ZNC}-TDC_{ZNA}",120,-30.,30.);
  fhTDCZNDiff->GetXaxis()->SetTitle("TDC_{ZNC}-TDC_{ZNA} (ns)");
  fOutput->Add(fhTDCZNDiff);     
  
  fhZNCSumQ = new TH1F("fhZNCSumQ", "ZNC sumPMQ_i signal", 250,0., 2000.);
  fOutput->Add(fhZNCSumQ);      
  fhZNASumQ = new TH1F("fhZNASumQ", "ZNA sumPMQ_i signal", 250,0., 2000.) ;
  fOutput->Add(fhZNASumQ);      
  fhZPCSumQ = new TH1F("fhZPCSumQ", "ZPC sumPMQ_i signal", 200,0., 2000.) ;
  fOutput->Add(fhZPCSumQ);      
  fhZPASumQ = new TH1F("fhZPASumQ", "ZPA sumPMQ_i signal", 200,0., 2000.) ;
  fOutput->Add(fhZPASumQ);      
  fhZEM1Spectrum = new TH1F("fhZEM1Spectrum", "ZEM1 signal", 100,0., 2500.);
  fOutput->Add(fhZEM1Spectrum);      
  fhZEM2Spectrum = new TH1F("fhZEM2Spectrum", "ZEM2 signal", 100,0., 2500.);
  fOutput->Add(fhZEM2Spectrum);      
  
  fhZNCpmc = new TH1F("fhZNCpmc","ZNC PMC", 200, 0., 2000.);
  fOutput->Add(fhZNCpmc);      
  fhZNApmc = new TH1F("fhZNApmc","ZNA PMC", 200, 0., 2000.); 
  fOutput->Add(fhZNApmc);      
  fhZPCpmc = new TH1F("fhZPCpmc","ZPC PMC", 200, 0., 2000.); 
  fOutput->Add(fhZPCpmc);      
  fhZPApmc = new TH1F("fhZPApmc","ZPA PMC", 200, 0., 2000.); 
  fOutput->Add(fhZPApmc);      
  
  fhZNCCentroid = new TH2F("fhZNCCentroid","Centroid over ZNC",70,-3.5,3.5,70,-3.5,3.5); 
  fOutput->Add(fhZNCCentroid);      
  fhZNACentroid = new TH2F("fhZNACentroid","Centroid over ZNA",70,-3.5,3.5,70,-3.5,3.5); 
  fOutput->Add(fhZNACentroid);      
  
  fhPMCZNCemd = new TH1F("fhPMCZNCemd","ZNC PMC lg",200, 0., 600.);   
  fOutput->Add(fhPMCZNCemd);      
  fhPMCZNAemd = new TH1F("fhPMCZNAemd","ZNA PMC lg",200, 0., 600.);   
  fOutput->Add(fhPMCZNAemd);     
  
  fDebunch = new TH2F("fDebunch","ZN TDC sum vs. diff", 120,-30,30,120,-30,-30);
  fOutput->Add(fDebunch);     
  
  fhTDCZNAcorr = new TH1F("fhTDCZNAcorr","ZNA TDC corrected",120, -20., 20.);
  fOutput->Add(fhTDCZNAcorr);
  
  fhTDCZNCcorr = new TH1F("fhTDCZNCcorr","ZNC TDC corrected",120, -20., 20.);
  fOutput->Add(fhTDCZNCcorr);
    
  
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskZDCpA::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  if(fDebug>1) printf(" **** AliAnalysisTaskZDCpA::UserExec() \n");
  
  if (!InputEvent()) {
    Printf("ERROR: InputEvent not available");
    return;
  }

      
  AliESDEvent* esd = dynamic_cast<AliESDEvent*> (InputEvent());
  if(!esd) return;
  // Select PHYSICS events (type=7, for data)
  //if(!fIsMCInput && esd->GetEventType()!=7) return; 
  
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

  }
  // ****************************************************
  
//  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    
  AliESDZDC *esdZDC = esd->GetESDZDC();
  
//  if((((AliInputEventHandler*)(am->GetInputEventHandler()))->IsEventSelected())){

    fhZEM1Spectrum->Fill(esdZDC->GetZDCEMEnergy(0));
    fhZEM2Spectrum->Fill(esdZDC->GetZDCEMEnergy(1));
  
    const Double_t * towZNC = esdZDC->GetZN1TowerEnergy();
    const Double_t * towZPC = esdZDC->GetZP1TowerEnergy();
    const Double_t * towZNA = esdZDC->GetZN2TowerEnergy();
    const Double_t * towZPA = esdZDC->GetZP2TowerEnergy();
    //    
    fhZNCpmc->Fill(towZNC[0]);    
    fhZNApmc->Fill(towZNA[0]);    
    fhZPCpmc->Fill(towZPC[0]);    
    fhZPApmc->Fill(towZPA[0]);
    //
    Float_t sumQZNC=0., sumQZNA=0., sumQZPC=0., sumQZPA=0.;
    for(int iq=1; iq<5; iq++){
       sumQZNC += towZNC[iq];
       sumQZNA += towZNA[iq];
       sumQZPC += towZPC [iq];
       sumQZPA += towZNC[iq];
    }     
  
    fhZNCSumQ->Fill(sumQZNC);	  
    fhZNASumQ->Fill(sumQZNA);
    fhZPCSumQ->Fill(sumQZPC);		  
    fhZPASumQ->Fill(sumQZPA);	  
  
    Double_t xyZNC[2]={-99.,-99.}, xyZNA[2]={-99.,-99.};
    esdZDC->GetZNCentroidInpp(xyZNC, xyZNA);
    
    fhZNCCentroid->Fill(xyZNC[0], xyZNC[1]); 
    fhZNACentroid->Fill(xyZNA[0], xyZNA[1]); 
    
    const Double_t * towZNCLG = esdZDC->GetZN1TowerEnergyLR();
    const Double_t * towZNALG = esdZDC->GetZN2TowerEnergyLR();
    fhPMCZNCemd->Fill(towZNCLG[0]);   
    fhPMCZNAemd->Fill(towZNALG[0]);   
  
//  }
    
  for(int i=0; i<4; i++){
    Float_t tdcC = 0.025*(esdZDC->GetZDCTDCData(10,i));
    Float_t tdcA = 0.025*(esdZDC->GetZDCTDCData(12,i));
    Float_t tdcL0 = 0.025*(esdZDC->GetZDCTDCData(15,i));
    if(tdcC != 0.) fhTDCZNCcorr->Fill(esdZDC->GetZDCTDCCorrected(10,i));
    if(tdcA != 0.) fhTDCZNAcorr->Fill(esdZDC->GetZDCTDCCorrected(12,i));
    //
    if(tdcC != 0.){
      fhTDCZNC->Fill(tdcC-tdcL0);
      if(tdcA != 0.){
        fhTDCZNA->Fill(tdcA-tdcL0);
        fhTDCZNDiff->Fill(esdZDC->GetZDCTDCCorrected(10,i)-esdZDC->GetZDCTDCCorrected(12,i));
        fhTDCZNSum->Fill(esdZDC->GetZDCTDCCorrected(10,i)-esdZDC->GetZDCTDCCorrected(12,i));
	fDebunch->Fill(esdZDC->GetZDCTDCCorrected(10,i)-esdZDC->GetZDCTDCCorrected(12,i),
	               esdZDC->GetZDCTDCCorrected(10,i)+esdZDC->GetZDCTDCCorrected(12,i));
      }
    }
  }  
  
  PostData(1, fOutput);
   
}



//________________________________________________________________________
void AliAnalysisTaskZDCpA::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
}
