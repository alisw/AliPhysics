#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskGammaJet.h"

#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDInputHandler.h"



#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODCaloCluster.h"
#include "AliGammaConversionAODObject.h"

#include <iostream>

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing
// Reviewed: A.Gheata (19/02/10)

using namespace std;

ClassImp(AliAnalysisTaskGammaJet)

//________________________________________________________________________
AliAnalysisTaskGammaJet::AliAnalysisTaskGammaJet(const char *name) : AliAnalysisTaskSE(name), 
  fESD(0), 
  fOutputList(0), 
  fHistPt(0),
  fHistPtPhos(0),
  fHistPtEmcal(0),
  fHistPtJets(0),
  fDeltaAODFileName("")
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskGammaJet::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputList = new TList();

  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 150, 0.1, 50);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
  
  fHistPtPhos = new TH1F("fHistPtPhos", "P_{T} distribution", 150, 0.1, 50);
  fHistPtPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtPhos->SetMarkerStyle(kFullCircle);
  
  fHistPtEmcal = new TH1F("fHistPtEmcal", "P_{T} distribution", 150, 0.1, 50);
  fHistPtEmcal->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtEmcal->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtEmcal->SetMarkerStyle(kFullCircle);


  fHistPtJets = new TH1F("fHistPtJets", "P_{T} distribution", 150, 0.1, 50);
  fHistPtJets->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtJets->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtJets->SetMarkerStyle(kFullCircle);
  
  fOutputList->Add(fHistPt);
  fOutputList->Add(fHistPtPhos);
  fOutputList->Add(fHistPtEmcal);
  fOutputList->Add(fHistPtJets);
}

//________________________________________________________________________
void AliAnalysisTaskGammaJet::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event



  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }

  TClonesArray * convGamma = NULL;
  AliAODEvent * aodEvent = AODEvent();
  
  if(!aodEvent) return;

  //Get Conversion gamma branch of AOD. 
  if(fDeltaAODFileName.Length() > 0) 
    convGamma = GetConversionGammas();
  else
    convGamma = dynamic_cast<TClonesArray*>(aodEvent->FindListObject("GammaConv"));

  if(convGamma) {
    for (Int_t iPhot = 0; iPhot < convGamma->GetEntriesFast(); iPhot++) {
      AliGammaConversionAODObject * aodO = dynamic_cast<AliGammaConversionAODObject*>(convGamma->At(iPhot));
      if (!aodO) {
	//printf("ERROR: Could not receive ga %d\n", iPhot);
	continue;
      }
      fHistPt->Fill(aodO->E());
    } 
  }

  if(fESD) {
    for(int icl = 0; icl < fESD->GetNumberOfCaloClusters(); icl++) {
      AliESDCaloCluster * cluster = fESD->GetCaloCluster(icl);
      if(cluster->GetEmcCpvDistance() > 10) {
	if (cluster->IsPHOS()) fHistPtPhos->Fill(cluster->E());
	if (cluster->IsEMCAL()) fHistPtEmcal->Fill(cluster->E());
      }
    }
  }

  if(aodEvent) {
    TClonesArray * jets = aodEvent->GetJets();
    if (jets) {
      for(int ij = 0; ij < jets->GetEntriesFast(); ij++) {
	AliAODJet * jet = aodEvent->GetJet(ij);
	if(jet) {
	  fHistPtJets->Fill(jet->E());
	}
      }
    }
  }
    
  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskGammaJet::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}

//________________________________________________________________________
TClonesArray * AliAnalysisTaskGammaJet::GetConversionGammas() {

  if(  !(fDeltaAODFileName.Length() > 0)  ) return NULL;
  
  //If AOD not in standard file have to locate it in extension
  AliAODHandler * aodHandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()); 
  if(aodHandler) {
    AliAODExtension * gExt = dynamic_cast<AliAODExtension*>(aodHandler->GetExtensions()->FindObject(fDeltaAODFileName));
    if(gExt) {
      AliAODEvent * aodEvent = gExt->GetAOD();
      return dynamic_cast<TClonesArray*>(aodEvent->FindListObject("GammaConv"));
    }
  }  
  return NULL;
}
