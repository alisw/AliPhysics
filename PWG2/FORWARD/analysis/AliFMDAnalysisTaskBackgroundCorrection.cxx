 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>
#include "TH2F.h"
#include "AliFMDAnalysisTaskBackgroundCorrection.h"
#include "AliAnalysisManager.h"
#include "AliESDFMD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliESDVertex.h"
#include "TMath.h"
#include "AliFMDAnaParameters.h"
//#include "AliFMDGeometry.h"

ClassImp(AliFMDAnalysisTaskBackgroundCorrection)


AliFMDAnalysisTaskBackgroundCorrection::AliFMDAnalysisTaskBackgroundCorrection()
: fDebug(0),
  fOutputList(0),
  fInputList(0),
  fHitList(0),
  fVertexString(0x0),
  fNevents(),
  fStandalone(kTRUE),
  fOutputVertexString(0)
{
  // Default constructor
  DefineInput (0, TList::Class());
  DefineOutput(0, TList::Class());

}
//_____________________________________________________________________
AliFMDAnalysisTaskBackgroundCorrection::AliFMDAnalysisTaskBackgroundCorrection(const char* name, Bool_t SE):
    AliAnalysisTask(name, "Density"),
    fDebug(0),
    fOutputList(0),
    fInputList(0),
    fHitList(0),
    fVertexString(0x0),
    fNevents(),
    fStandalone(kTRUE),
    fOutputVertexString(0)
{
  fStandalone = SE;
  if(fStandalone) {
    DefineInput (0, TList::Class());
    DefineOutput(0, TList::Class());
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBackgroundCorrection::CreateOutputObjects()
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  if(!fOutputList)
    fOutputList = new TList();
  fOutputList->SetName("BackgroundCorrectedPerEvent");
  if(!fHitList)
    fHitList = new TList();
  fHitList->SetName("HitsList");
  
  //if(fStandalone) {
  fOutputVertexString = new TObjString();
  // }
  fOutputList->Add(fOutputVertexString);
  
  
  
  TH2F* hMult = 0;
  TH2F* hMultTrVtx = 0;
  TH2F* hHits = 0;
  // TH2F* hHitsNoCuts = 0;
  Int_t nVtxbins = pars->GetNvtxBins();
  
  for(Int_t det =1; det<=3;det++)
    {
      Int_t nRings = (det==1 ? 1 : 2);
      for(Int_t ring = 0;ring<nRings;ring++)
	{
	  Char_t ringChar = (ring == 0 ? 'I' : 'O');
	  Int_t  nSec     = (ring == 0 ? 20 : 40);
	  
	  for(Int_t i = 0; i< nVtxbins; i++) {
	    TH2F* hBg = pars->GetBackgroundCorrection(det, ringChar, i);
	    hMult  = new TH2F(Form("mult_FMD%d%c_vtxbin%d",det,ringChar,i),Form("mult_FMD%d%c_vtxbin%d",det,ringChar,i),
			      hBg->GetNbinsX(),
			      hBg->GetXaxis()->GetXmin(),
			      hBg->GetXaxis()->GetXmax(),
			      nSec, 0, 2*TMath::Pi());
	    hMult->Sumw2();
	    fOutputList->Add(hMult);
	    hMultTrVtx  = new TH2F(Form("multTrVtx_FMD%d%c_vtxbin%d",det,ringChar,i),Form("mult_FMD%d%c_vtxbin%d",det,ringChar,i),
				   hBg->GetNbinsX(),
				   hBg->GetXaxis()->GetXmin(),
				   hBg->GetXaxis()->GetXmax(),
				   nSec, 0, 2*TMath::Pi());
	    hMultTrVtx->Sumw2();

	    fOutputList->Add(hMultTrVtx);
	    hHits  = new TH2F(Form("hits_FMD%d%c_vtxbin%d",det,ringChar,i),Form("hits_FMD%d%c_vtxbin%d",det,ringChar,i),
			      hBg->GetNbinsX(),
			      hBg->GetXaxis()->GetXmin(),
			      hBg->GetXaxis()->GetXmax(),
			      nSec, 0, 2*TMath::Pi());
	    
	    /*  hHitsNoCuts  = new TH2F(Form("hits_NoCuts_FMD%d%c_vtxbin%d",det,ringChar,i),Form("hits_NoCuts_FMD%d%c_vtxbin%d",det,ringChar,i),
				    hBg->GetNbinsX(),
				    hBg->GetXaxis()->GetXmin(),
				    hBg->GetXaxis()->GetXmax(),
				    nSec, 0, 2*TMath::Pi());
	    
	    */
	    hHits->Sumw2();
	    fHitList->Add(hHits);
	    fOutputList->Add(hHits);
	    	    
	  }
	} 
    }
  
  
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBackgroundCorrection::ConnectInputData(Option_t */*option*/)
{
  if(fStandalone) {
    fInputList   = (TList*)GetInputData(0);
    
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBackgroundCorrection::Exec(Option_t */*option*/)
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  fVertexString = (TObjString*)fInputList->At(0);
   
  Int_t vtxbin   = fVertexString->GetString().Atoi();
  fOutputVertexString->SetString(Form("%d",vtxbin));
  
  fNevents.Fill(vtxbin);
  //Reset everything
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
      TH2F* hMult = (TH2F*)fOutputList->FindObject(Form("mult_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      hMult->Reset();
      TH2F* hMultTrVtx = (TH2F*)fOutputList->FindObject(Form("multTrVtx_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      hMultTrVtx->Reset();
    }
    
  }
  
  
  
  for(UShort_t det=1;det<=3;det++) {
    
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
      
      TH2F* hMult      = (TH2F*)fOutputList->FindObject(Form("mult_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hMultTrVtx = (TH2F*)fOutputList->FindObject(Form("multTrVtx_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hMultInput = (TH2F*)fInputList->FindObject(Form("FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hHits      = (TH2F*)fOutputList->FindObject(Form("hits_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      
      if(pars->GetProcessHits())
	hHits->Add(hMultInput);
      
      TH2F* hBg        = pars->GetBackgroundCorrection(det, ringChar, vtxbin);
      
      hMult->Add(hMultInput);
      hMultTrVtx->Add(hMultInput);
      hMult->Divide(hBg);//,"B");
      hMultTrVtx->Divide(hBg);//,"B");

      //sharing efficiency correction ?
      
      TH1F* hSharingEff = pars->GetSharingEfficiency(det,ringChar,vtxbin);
      TH1F* hSharingEffTrVtx = pars->GetSharingEfficiencyTrVtx(det,ringChar,vtxbin);	
      
      for(Int_t nx=1; nx<hMult->GetNbinsX(); nx++) {
	Float_t correction = hSharingEff->GetBinContent(nx);
	Float_t correctionTrVtx = hSharingEffTrVtx->GetBinContent(nx);
	for(Int_t ny=1; ny<hMult->GetNbinsY(); ny++) {
	  
	  if(correction != 0){
	    hMult->SetBinContent(nx,ny,hMult->GetBinContent(nx,ny)/correction);
	    Float_t error = TMath::Sqrt(TMath::Power(hMult->GetBinError(nx,ny),2) + TMath::Power(hMult->GetBinContent(nx,ny)*hSharingEff->GetBinError(nx),2)) / correction;
	    hMult->SetBinError(nx,ny,error);
	  }
	  if(correctionTrVtx != 0){
	    hMultTrVtx->SetBinContent(nx,ny,hMultTrVtx->GetBinContent(nx,ny)/correctionTrVtx);
	    Float_t error = TMath::Sqrt(TMath::Power(hMultTrVtx->GetBinError(nx,ny),2) + TMath::Power(hMultTrVtx->GetBinContent(nx,ny)*hSharingEffTrVtx->GetBinError(nx),2)) / correctionTrVtx;
	    hMultTrVtx->SetBinError(nx,ny,error);
	  }
	}
	
      }
      
      hMult->Scale(1/pars->GetEventSelectionEfficiency(vtxbin));
      
      
    }
  }
  if(fStandalone) {
    PostData(0, fOutputList); 
  }
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBackgroundCorrection::Terminate(Option_t */*option*/) {
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  Int_t nVtxbins = pars->GetNvtxBins();
  
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
      for(Int_t i =0; i<nVtxbins; i++) {
	TH2F* hHits      = (TH2F*)fOutputList->FindObject(Form("hits_FMD%d%c_vtxbin%d",det,ringChar,i));
	TH1D* hHitsproj  = hHits->ProjectionX(Form("hits_FMD%d%c_vtxbin%d_proj",det,ringChar,i),1,hHits->GetNbinsY());
	TH1D* hHitsNoCuts = (TH1D*)hHitsproj->Clone(Form("hits_NoCuts_FMD%d%c_vtxbin%d_proj",det,ringChar,i));
	
	hHitsNoCuts->Scale(1/pars->GetEventSelectionEfficiency(i));
	fHitList->Add(hHitsproj);
	fHitList->Add(hHitsNoCuts);
	
      }
    }
  }
}
//_____________________________________________________________________
//
//
// EOF
