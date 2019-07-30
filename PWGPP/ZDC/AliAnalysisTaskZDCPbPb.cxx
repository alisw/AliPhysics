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
#include "AliAnalysisTaskZDCPbPb.h"

ClassImp(AliAnalysisTaskZDCPbPb)


//________________________________________________________________________
AliAnalysisTaskZDCPbPb::AliAnalysisTaskZDCPbPb():
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
    fhZNCpmcUncalib(0x0),
    fhZNApmcUncalib(0x0),
    fhZPCpmcUncalib(0x0),
    fhZPApmcUncalib(0x0),
    fhZNCpmc(0x0),
    fhZNApmc(0x0),
    fhZPCpmc(0x0),
    fhZPApmc(0x0),
    fhZNCCentroid(0x0),
    fhZNACentroid(0x0),
    fhPMCZNCemdUncalib(0x0),
    fhPMCZNAemdUncalib(0x0),
    fhPMCZNCemd(0x0),
    fhPMCZNAemd(0x0),
    fDebunch(0x0),
    fhTDCZNC(0x0),
    fhTDCZNA(0x0)
{
   // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskZDCPbPb::AliAnalysisTaskZDCPbPb(const char *name):
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
    fhZNCpmcUncalib(0x0),
    fhZNApmcUncalib(0x0),
    fhZPCpmcUncalib(0x0),
    fhZPApmcUncalib(0x0),
    fhZNCpmc(0x0),
    fhZNApmc(0x0),
    fhZPCpmc(0x0),
    fhZPApmc(0x0),
    fhZNCCentroid(0x0),
    fhZNACentroid(0x0),
    fhPMCZNCemdUncalib(0x0),
    fhPMCZNAemdUncalib(0x0),
    fhPMCZNCemd(0x0),
    fhPMCZNAemd(0x0),
    fDebunch(0x0) ,
    fhTDCZNC(0x0),
    fhTDCZNA(0x0)
{
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());

}

//________________________________________________________________________
AliAnalysisTaskZDCPbPb& AliAnalysisTaskZDCPbPb::operator=(const AliAnalysisTaskZDCPbPb& c)
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
AliAnalysisTaskZDCPbPb::AliAnalysisTaskZDCPbPb(const AliAnalysisTaskZDCPbPb& ana):
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
  fhZNCpmcUncalib(ana.fhZNCpmcUncalib),
  fhZNApmcUncalib(ana.fhZNApmcUncalib),
  fhZPCpmcUncalib(ana.fhZPCpmcUncalib),
  fhZPApmcUncalib(ana.fhZPApmcUncalib),
  fhZNCpmc(ana.fhZNCpmc),
  fhZNApmc(ana.fhZNApmc),
  fhZPCpmc(ana.fhZPCpmc),
  fhZPApmc(ana.fhZPApmc),
  fhZNCCentroid(ana.fhZNCCentroid),
  fhZNACentroid(ana.fhZNACentroid),
  fhPMCZNCemdUncalib(ana.fhPMCZNCemdUncalib),
  fhPMCZNAemdUncalib(ana.fhPMCZNAemdUncalib),
  fhPMCZNCemd(ana.fhPMCZNCemd),
  fhPMCZNAemd(ana.fhPMCZNAemd),
  fDebunch(ana.fDebunch),
  fhTDCZNC(ana.fhTDCZNC),
  fhTDCZNA(ana.fhTDCZNA)
{
  //
  // Copy Constructor
  //
}

//________________________________________________________________________
AliAnalysisTaskZDCPbPb::~AliAnalysisTaskZDCPbPb()
{
  // Destructor
  if(fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    delete fOutput; fOutput=0;
  }

}

//________________________________________________________________________
void AliAnalysisTaskZDCPbPb::UserCreateOutputObjects()
{
  // Create the output containers

  fOutput = new TList;
  fOutput->SetOwner();
  //fOutput->SetName("output");

  // Setting different upper edges for MC and for data
  Float_t upUncalZN = 1400., upCalZN = 250.;
  Float_t upUncalZP = 1200., upCalZP = 100.;
  Float_t upZEM = 1400.;
  if(fIsMCInput){
    upUncalZN = 4000.;
    upCalZN = 100;
    upUncalZP = 4000.;
    upCalZP = 80;
    upZEM = 800.;
  }

  fDebunch = new TH2F("fDebunch","ZN TDC sum vs. diff", 120,-30,30,120,-30,30);
  fOutput->Add(fDebunch);

  fhTDCZNSum = new TH1F("fhTDCZNSum","TDC_{ZNC}+TDC_{ZNA}",120,-30.,30.);
  fhTDCZNSum->GetXaxis()->SetTitle("TDC_{ZNC}+TDC_{ZNA} (ns)");
  fOutput->Add(fhTDCZNSum);

  fhTDCZNDiff = new TH1F("fhTDCZNDiff","TDC_{ZNC}-TDC_{ZNA}",120,-30.,30.);
  fhTDCZNDiff->GetXaxis()->SetTitle("TDC_{ZNC}-TDC_{ZNA} (ns)");
  fOutput->Add(fhTDCZNDiff);

  fhZNCSpectrum = new TH1F("fhZNCSpectrum", "ZNC signal", 250,0., upCalZN);
  fOutput->Add(fhZNCSpectrum);
  fhZNASpectrum = new TH1F("fhZNASpectrum", "ZNA signal", 250,0., upCalZN) ;
  fOutput->Add(fhZNASpectrum);
  fhZPCSpectrum = new TH1F("fhZPCSpectrum", "ZPC signal", 150,0., upCalZP) ;
  fOutput->Add(fhZPCSpectrum);
  fhZPASpectrum = new TH1F("fhZPASpectrum", "ZPA signal", 150,0., upCalZP) ;
  fOutput->Add(fhZPASpectrum);
  fhZEM1Spectrum = new TH1F("fhZEM1Spectrum", "ZEM1 signal", 200,0., upZEM);
  fOutput->Add(fhZEM1Spectrum);
  fhZEM2Spectrum = new TH1F("fhZEM2Spectrum", "ZEM2 signal", 200,0., upZEM);
  fOutput->Add(fhZEM2Spectrum);

  fhZNCpmcUncalib = new TH1F("fhZNCpmcUncalib","ZNC PMC NO ENERGY calibration",200, 0., upUncalZN);
  fhZNCpmcUncalib->SetXTitle("ZNC signa (ADC channels)");
  fOutput->Add(fhZNCpmcUncalib);
  fhZNApmcUncalib = new TH1F("fhZNApmcUncalib","ZNA PMC NO ENERGY calibration",200, 0., upUncalZN);
  fhZNApmcUncalib->SetXTitle("ZNA signa (ADC channels)");
  fOutput->Add(fhZNApmcUncalib);
  fhZPCpmcUncalib = new TH1F("fhZPCpmcUncalib","ZPC PMC NO ENERGY calibration",200, 0., upUncalZP);
  fhZPCpmcUncalib->SetXTitle("ZPC signa (ADC channels)");
  fOutput->Add(fhZPCpmcUncalib);
  fhZPApmcUncalib = new TH1F("fhZPApmcUncalib","ZPA PMC NO ENERGY calibration",200, 0., upUncalZP);
  fhZPApmcUncalib->SetXTitle("ZPA signa (ADC channels)");
  fOutput->Add(fhZPApmcUncalib);

  fhZNCpmc = new TH1F("fhZNCpmc","ZNC PMC",250, 0., upCalZN);
  fhZNCpmc->SetXTitle("ZNC energy (TeV)");
  fOutput->Add(fhZNCpmc);
  fhZNApmc = new TH1F("fhZNApmc","ZNA PMC",250, 0., upCalZN);
  fhZNApmc->SetXTitle("ZNA energy (TeV)");
  fOutput->Add(fhZNApmc);
  fhZPCpmc = new TH1F("fhZPCpmc","ZPC PMC",100, 0., upCalZP);
  fhZPCpmc->SetXTitle("ZPC energy (TeV)");
  fOutput->Add(fhZPCpmc);
  fhZPApmc = new TH1F("fhZPApmc","ZPA PMC",100, 0., upCalZP);
  fhZPApmc->SetXTitle("ZPA energy (TeV)");
  fOutput->Add(fhZPApmc);

  fhZNCCentroid = new TH2F("fhZNCCentroid","Centroid over ZNC",70,-3.5,3.5,70,-3.5,3.5);
  fOutput->Add(fhZNCCentroid);
  fhZNACentroid = new TH2F("fhZNACentroid","Centroid over ZNA",70,-3.5,3.5,70,-3.5,3.5);
  fOutput->Add(fhZNACentroid);

  fhPMCZNCemdUncalib = new TH1F("fhPMCZNCemdUncalib","ZNC PMC lg",200, 0., 10.);
  fOutput->Add(fhPMCZNCemdUncalib);
  fhPMCZNAemdUncalib = new TH1F("fhPMCZNAemdUncalib","ZNA PMC lg",200, 0., 10.);
  fOutput->Add(fhPMCZNAemdUncalib);

  fhPMCZNCemd = new TH1F("fhPMCZNCemd","ZNC PMC lg",200, 0., 800.);
  fOutput->Add(fhPMCZNCemd);
  fhPMCZNAemd = new TH1F("fhPMCZNAemd","ZNA PMC lg",200, 0., 800.);
  fOutput->Add(fhPMCZNAemd);

  fhTDCZNC = new TH1F("fhTDCZNC","TDC_{ZNC}",60,-30.,30.);
  fhTDCZNC->GetXaxis()->SetTitle("TDC_{ZNC} (ns)");
  fOutput->Add(fhTDCZNC);

  fhTDCZNA = new TH1F("fhTDCZNA","TDC_{ZNA}",60,-30.,30.);
  fhTDCZNA->GetXaxis()->SetTitle("TDC_{ZNA} (ns)");
  fOutput->Add(fhTDCZNA);

  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskZDCPbPb::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  if(fDebug>1) printf(" **** AliAnalysisTaskZDCPbPb::UserExec() \n");

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

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();

  AliESDZDC *esdZDC = esd->GetESDZDC();

  if((((AliInputEventHandler*)(am->GetInputEventHandler()))->IsEventSelected())){

    fhZNCSpectrum->Fill(esdZDC->GetZDCN1Energy());
    fhZNASpectrum->Fill(esdZDC->GetZDCN2Energy());
    fhZPCSpectrum->Fill(esdZDC->GetZDCP1Energy());
    fhZPASpectrum->Fill(esdZDC->GetZDCP2Energy());
    fhZEM1Spectrum->Fill(esdZDC->GetZDCEMEnergy(0));
    fhZEM2Spectrum->Fill(esdZDC->GetZDCEMEnergy(1));

    const Double_t * towZNC = esdZDC->GetZN1TowerEnergy();
    const Double_t * towZPC = esdZDC->GetZP1TowerEnergy();
    const Double_t * towZNA = esdZDC->GetZN2TowerEnergy();
    const Double_t * towZPA = esdZDC->GetZP2TowerEnergy();
    //
    fhZNCpmcUncalib->Fill(towZNC[0]);
    fhZNApmcUncalib->Fill(towZNA[0]);
    fhZPCpmcUncalib->Fill(towZPC[0]);
    fhZPApmcUncalib->Fill(towZPA[0]);
    //
    fhZNCpmc->Fill(towZNC[0]/1000.);
    fhZNApmc->Fill(towZNA[0]/1000.);
    fhZPCpmc->Fill(towZPC[0]/1000.);
    fhZPApmc->Fill(towZPA[0]/1000.);

    Double_t xyZNC[2]={-99.,-99.}, xyZNA[2]={-99.,-99.};
    esdZDC->GetZNCentroidInPbPb(2510., xyZNC, xyZNA);
    //esdZDC->GetZNCentroidInpp(xyZNC, xyZNA);

    fhZNCCentroid->Fill(xyZNC[0], xyZNC[1]);
    fhZNACentroid->Fill(xyZNA[0], xyZNA[1]);
  }

  const Double_t * towZNCLG = esdZDC->GetZN1TowerEnergyLR();
  const Double_t * towZNALG = esdZDC->GetZN2TowerEnergyLR();
  fhPMCZNCemdUncalib->Fill(towZNCLG[0]);
  fhPMCZNAemdUncalib->Fill(towZNALG[0]);

  fhPMCZNCemd->Fill(towZNCLG[0]/1000.);
  fhPMCZNAemd->Fill(towZNALG[0]/1000.);

  Float_t tdcC=999., tdcA=999;
  Float_t tdcSum=999., tdcDiff=999;
  for(int i=0; i<4; i++){
    if(esdZDC->GetZDCTDCData(esdZDC->GetZNCTDCChannel(),i) != 0.){
      tdcC = esdZDC->GetZDCTDCCorrected(esdZDC->GetZNCTDCChannel(),i);
      fhTDCZNC->Fill(tdcC);
      if(esdZDC->GetZDCTDCData(esdZDC->GetZNATDCChannel(),i) != 0.){
        tdcA = esdZDC->GetZDCTDCCorrected(esdZDC->GetZNATDCChannel(),i);
        fhTDCZNA->Fill(tdcA);
        tdcSum = tdcC+tdcA;
        tdcDiff = tdcC-tdcA;
      }
    }
    if(tdcSum!=999.) fhTDCZNSum->Fill(tdcSum);
    if(tdcDiff!=999.)fhTDCZNDiff->Fill(tdcDiff);
    if(tdcSum!=999. && tdcDiff!=999.)  fDebunch->Fill(tdcDiff, tdcSum);
  }

  PostData(1, fOutput);

}



//________________________________________________________________________
void AliAnalysisTaskZDCPbPb::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
}
