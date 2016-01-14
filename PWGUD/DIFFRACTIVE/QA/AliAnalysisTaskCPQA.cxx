// QA task for central production study
// author: Martin Poghosyan

#include <TList.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>


#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliHeader.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultiplicity.h"
#include "AliTriggerAnalysis.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisTaskCPQA.h"


ClassImp(AliAnalysisTaskCPQA)

//________________________________________________________________________
  AliAnalysisTaskCPQA::AliAnalysisTaskCPQA(const char *name)
  : AliAnalysisTaskSE(name),
  fUseMC(kFALSE),
  fESD(0),
  fOutputList(0),
  fhEvent(0),
  fTriggerAnalysis(0)
{
  for(Int_t i = 0; i<4; i++)
    {
      fhV0A[i]=0;
      fhV0C[i]=0;
      fhV0online[i]=0;
      fhV0offline[i]=0;
      fhSPDFiredChip[i]=0;
      fhSPDFastOrChip[i]=0;
      fhReferenceMultiplicity[i]=0;
      fhVtxTrack[i]=0;
    }

// fEtaMaxM = 2;
// fEtaMaxD = 0.9;
// fVtxZmax = 10;

  DefineOutput(1, TList::Class());
}

AliAnalysisTaskCPQA::~AliAnalysisTaskCPQA() 
{
   if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    printf("Deleteing output\n");

    if(fOutputList){
      delete fOutputList;
      fOutputList = 0;
    }

        if(fTriggerAnalysis)
          delete fTriggerAnalysis;

   }
}

//________________________________________________________________________
void AliAnalysisTaskCPQA::UserCreateOutputObjects()
{
  fTriggerAnalysis = new AliTriggerAnalysis();
  if (fUseMC) fTriggerAnalysis->SetAnalyzeMC(1);

  // Create histograms
  // Called once

  fOutputList = new TList();
  fOutputList->SetOwner();


  fhEvent = new TH1F("hEvent","hEvent",100, -0.5, 99.5); 
  fOutputList->Add(fhEvent);

  for(Int_t i = 0; i<4; i++)
    {
      fhV0A[i] = Hist2D(Form("hV0A_%d",i), 5, -1.5, 3.5, 5, -1.5, 3.5,"V0A_{online}","V0A_{offline}"); fOutputList->Add(fhV0A[i]);
      fhV0C[i] = Hist2D(Form("hV0C_%d",i), 5, -1.5, 3.5, 5, -1.5, 3.5,"V0C_{online}","V0C_{offline}"); fOutputList->Add(fhV0C[i]);
      fhV0online[i]  = Hist2D(Form("hV0online_%d",i) , 5, -1.5, 3.5, 5, -1.5, 3.5,"V0C_{online}","V0A_{online}"); fOutputList->Add(fhV0online[i]);
      fhV0offline[i] = Hist2D(Form("hV0offline_%d",i), 5, -1.5, 3.5, 5, -1.5, 3.5,"V0C_{offline}","V0A_{offline}"); fOutputList->Add(fhV0offline[i]);
      fhSPDFiredChip[i] = Hist1D(Form("fhSPDFiredChip_%d",i), 1200, -0.5, 1199.5); fOutputList->Add(fhSPDFiredChip[i]);
      fhSPDFastOrChip[i] = Hist1D(Form("fhSPDFastOrChip_%d",i), 1200, -0.5, 1199.5); fOutputList->Add(fhSPDFastOrChip[i]);
      fhReferenceMultiplicity[i] = Hist1D(Form("fhReferenceMultiplicity_%d",i), 50, -10.5, 39.5);  fOutputList->Add(fhReferenceMultiplicity[i]);
      fhVtxTrack[i] = Hist3D(Form("fhVtxTrack_%d",i), 100, -1, 1, 100, -1, 1, 1000, -30, 30, "x_{vtx}", "y_{vtx}", "z_{vtx}"); fOutputList->Add(fhVtxTrack[i]);
   }

  PostData(1, fOutputList);
 

}



TH1F* AliAnalysisTaskCPQA::Hist1D(const char* name, Int_t nBins, Double_t xMin, Double_t xMax,  const char* xLabel, Int_t color, Int_t lst, const char* yLabel)
{
// create a 1D histogram

  TH1F* res = new TH1F(name, name, nBins, xMin, xMax);
  if (xLabel) res->GetXaxis()->SetTitle(xLabel);
  if (yLabel) res->GetYaxis()->SetTitle(yLabel);
  res->SetLineColor(color);
  res->SetMarkerColor(color);
  res->SetLineStyle(lst);
  return res;
}


TH2F *AliAnalysisTaskCPQA::Hist2D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t yMin, Double_t yMax, const char* xLabel, const char* yLabel, Int_t color)
{
// create a 2D histogram

  TH2F *res = new TH2F(name, name, nBinsx, xMin, xMax, nBinsy, yMin, yMax);
  if (xLabel) res->GetXaxis()->SetTitle(xLabel);
  if (yLabel) res->GetYaxis()->SetTitle(yLabel);
  //  res->SetMarkerStyle(kFullCircle);
  //  res->SetOption("E");
  res->SetLineColor(color);
  //  fOutputList->Add(res);
  return res;
}

  TH3F *AliAnalysisTaskCPQA::Hist3D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t yMin, Double_t yMax,  Int_t nBinsz, Double_t zMin, Double_t zMax, const char* xLabel, const char* yLabel, const char *zLabel)
{
// create a 3D histogram

  TH3F *res = new TH3F(name, name, nBinsx, xMin, xMax, nBinsy, yMin, yMax, nBinsz, zMin, zMax);
  if (xLabel) res->GetXaxis()->SetTitle(xLabel);
  if (yLabel) res->GetYaxis()->SetTitle(yLabel);
  if (zLabel) res->GetZaxis()->SetTitle(zLabel);
  //  res->SetMarkerStyle(kFullCircle);
  //  res->SetOption("E");
  //  res->SetLineColor(color);
  //  fOutputList->Add(res);
  return res;
}



//________________________________________________________________________
void AliAnalysisTaskCPQA::UserExec(Option_t *)
{
  // Main loop
  // Called for each event

  //    return;
  AliVEvent *event = InputEvent();
  if (!event) {
     Error("UserExec", "Could not retrieve event");
     return;
  }

  fESD = dynamic_cast<AliESDEvent*> (InputEvent());
  if (fESD) {
    LoopESD();
      if (fUseMC)
         LoopESDMC();
  }



  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskCPQA::Terminate(Option_t*)
{
 
}


void AliAnalysisTaskCPQA::LoopESDMC()
{
  // Main loop
  // Called for each event
  /*

  Int_t indexD1 = 1;
  Int_t indexD2 = 2;


  AliMCEvent *mcEvent = (AliMCEvent*) MCEvent();
  if (!mcEvent) {
    Error("LoopESDMC", "Could not retrieve MC event");
    return;
  }


  AliHeader* header = mcEvent->Header();
  if (!header)
    {
          AliDebug(AliLog::kError, "Header not available");
          return;
    }

  AliGenEventHeader* genHeader = header->GenEventHeader();
  if(!genHeader)
    {
      AliDebug(AliLog::kError, "GenHeader not available");
      return;
    }

*/  
 

}




//________________________________________________________________________
void AliAnalysisTaskCPQA::LoopESD()
{
   Int_t TrType = 0;
   Bool_t fkIsPhysSel = kFALSE;

  fhEvent->Fill(0);

  if(!fUseMC)
    {
      TrType = -1;
      if(fESD->IsTriggerClassFired("CINT1-B-NOPF-ALLNOTRD"))        TrType=0;
      else if(fESD->IsTriggerClassFired("CINT1-AC-NOPF-ALLNOTRD"))  TrType=1;
      else if(fESD->IsTriggerClassFired("CINT1-E-NOPF-ALLNOTRD"))   TrType=3;

      else if(fESD->IsTriggerClassFired("CINT1-B-NOPF-FASTNOTRD"))  TrType=0;
      else if(fESD->IsTriggerClassFired("CINT1-AC-NOPF-FASTNOTRD")) TrType=1;
      else if(fESD->IsTriggerClassFired("CINT1-E-NOPF-FASTNOTRD"))  TrType=3;

      else if(fESD->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL"))  TrType=0;
      else if(fESD->IsTriggerClassFired("CINT1A-ABCE-NOPF-ALL"))  TrType=1;
      else if(fESD->IsTriggerClassFired("CINT1C-ABCE-NOPF-ALL"))  TrType=2;
      else if(fESD->IsTriggerClassFired("CINT1-E-NOPF-ALL"))      TrType=3;

      UInt_t mask =  ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
      fkIsPhysSel=(mask & AliVEvent::kMB) ? 1 : 0; // check if minimum bias trigger class fired
 

      if(!fkIsPhysSel) return;
      printf("TrType = %d\n",TrType);
      if(TrType==-1) return;
    }

  fhEvent->Fill(1);


  Int_t V0Aonline  = fTriggerAnalysis->V0Trigger(fESD, AliTriggerAnalysis::kASide, kTRUE) ;
  Int_t V0Conline  = fTriggerAnalysis->V0Trigger(fESD, AliTriggerAnalysis::kCSide, kTRUE) ;
  Int_t V0Aoffline = fTriggerAnalysis->V0Trigger(fESD, AliTriggerAnalysis::kASide, kFALSE) ;
  Int_t V0Coffline = fTriggerAnalysis->V0Trigger(fESD, AliTriggerAnalysis::kCSide, kFALSE) ;


  fhV0A[TrType]->Fill(V0Aonline, V0Aoffline);
  fhV0C[TrType]->Fill(V0Conline, V0Coffline);
  fhV0online[TrType]->Fill(V0Conline, V0Aonline);
  fhV0offline[TrType]->Fill(V0Coffline, V0Aoffline);


  if(V0Aoffline!=0 || V0Coffline!=0) 
    {
      return;
    }
 
  fhEvent->Fill(2);

  const AliMultiplicity *mult = fESD->GetMultiplicity();
  for (Int_t i=0; i<1200; i++)
    {
      if(mult->TestFiredChipMap(i)) fhSPDFiredChip[TrType]->Fill(i);
      if(mult->TestFastOrFiredChips(i)) fhSPDFastOrChip[TrType]->Fill(i);
    }


  fhReferenceMultiplicity[TrType]->Fill(AliESDtrackCuts::GetReferenceMultiplicity(fESD,AliESDtrackCuts::kTrackletsITSTPC,1.2));


  const AliESDVertex *primaryTrackVtx = fESD->GetPrimaryVertexTracks();
  if (!primaryTrackVtx->GetStatus() && !primaryTrackVtx->GetStatus() )
    {
      return;
    }

  fhEvent->Fill(3);

  fhVtxTrack[TrType]->Fill(primaryTrackVtx->GetX(), primaryTrackVtx->GetY(), primaryTrackVtx->GetZ());



  return;
}   

