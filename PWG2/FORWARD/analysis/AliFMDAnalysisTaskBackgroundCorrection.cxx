 
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
  TH2F* hHits = 0;
  TH2F* hHitsNoCuts = 0;
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
	    //hHitsNoCuts->Sumw2();
	    
	    fHitList->Add(hHits);
	    fOutputList->Add(hHits);
	    // fHitList->Add(hHitsNoCuts);
	    //  fOutputList->Add(hHitsNoCuts);
	    
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
    }
    
  }
  
  
  
  for(UShort_t det=1;det<=3;det++) {
   
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
   
      TH2F* hMultTotal = (TH2F*)fOutputList->FindObject(Form("mult_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
    
      TH2F* hMultInput = (TH2F*)fInputList->FindObject(Form("FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hHits      = (TH2F*)fOutputList->FindObject(Form("hits_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      
      if(pars->GetProcessHits())
	 hHits->Add(hMultInput);
      
      TH2F* hBg        = pars->GetBackgroundCorrection(det, ringChar, vtxbin);
      
      hMultTotal->Add(hMultInput);
      
      hMultTotal->Divide(hBg);//,"B");
      /*for(Int_t i = 1; i<=hTmp->GetNbinsX();i++) {
	for(Int_t j = 1; j<=hTmp->GetNbinsY();j++) {
	  Float_t mult = hTmp->GetBinContent(i,j);
	  if(mult == 0) continue;
	  Float_t correction = hBg->GetBinContent(i,j);
	  
	  Float_t multcor = mult;
	  if(correction != 0)
	    multcor = multcor/correction;
	  else
	    std::cout<<"Warning! No correction for bin "<<i<<" , "<<j<<std::endl;
	  
	  hTmp->SetBinContent(i,j,multcor);
	}
	}*/
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
