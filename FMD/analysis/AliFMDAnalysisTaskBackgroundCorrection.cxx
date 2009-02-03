 
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
#include "AliFMDGeometry.h"

ClassImp(AliFMDAnalysisTaskBackgroundCorrection)


AliFMDAnalysisTaskBackgroundCorrection::AliFMDAnalysisTaskBackgroundCorrection()
: fDebug(0),
  fOutputList(0),
  fInputList(0),
  fArray(0),
  fInputArray(0),
  fVertexString(0x0),
  fNevents(),
  fStandalone(kTRUE),
  fOutputVertexString(0)
{
  // Default constructor
  DefineInput (0, TList::Class());
  DefineOutput(0, TList::Class());
  DefineOutput(1, TObjString::Class());
}
//_____________________________________________________________________
AliFMDAnalysisTaskBackgroundCorrection::AliFMDAnalysisTaskBackgroundCorrection(const char* name, Bool_t SE):
    AliAnalysisTask(name, "Density"),
    fDebug(0),
    fOutputList(0),
    fInputList(0),
    fArray(),
    fInputArray(0),
    fVertexString(0x0),
    fNevents(),
    fStandalone(kTRUE),
    fOutputVertexString(0)
{
  fStandalone = SE;
  if(fStandalone) {
    DefineInput (0, TList::Class());
    DefineOutput(0, TList::Class());
    DefineOutput(1, TObjString::Class());
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBackgroundCorrection::CreateOutputObjects()
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  fArray.SetName("FMD");
  fArray.SetOwner();
  
  if(!fOutputList)
    fOutputList = new TList();
  fOutputList->SetName("BackgroundCorrectedPerEvent");
  
  
  TH2F* hMult = 0;
  TH2F* hHits = 0;
  Int_t nVtxbins = pars->GetNvtxBins();
  
  for(Int_t det =1; det<=3;det++)
    {
      TObjArray* detArray = new TObjArray();
      detArray->SetName(Form("FMD%d",det));
      fArray.AddAtAndExpand(detArray,det);
      Int_t nRings = (det==1 ? 1 : 2);
      for(Int_t ring = 0;ring<nRings;ring++)
	{
	  Char_t ringChar = (ring == 0 ? 'I' : 'O');
	  Int_t  nSec     = (ring == 0 ? 20 : 40);
	  
	  TObjArray* vtxArray = new TObjArray();
	  vtxArray->SetName(Form("FMD%d%c",det,ringChar));
	  detArray->AddAtAndExpand(vtxArray,ring);
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
	    
	    hHits->Sumw2();
	    fOutputList->Add(hHits);
	    vtxArray->AddAtAndExpand(hMult,i);
	    
	  }
	} 
    }
  if(fStandalone) {
    fOutputVertexString = new TObjString();
  }
  fOutputList->Add(fOutputVertexString);
  
  
  
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
  
  fInputArray   = (TObjArray*)fInputList->At(0);
  fVertexString = (TObjString*)fInputList->At(1);
  
  
  Int_t vtxbin   = fVertexString->GetString().Atoi();
  fOutputVertexString->SetString(Form("%d",vtxbin));
  
  //fNevents.operator[](vtxbin)++;
  fNevents.Fill(vtxbin);
  
  //Reset everything
  for(UShort_t det=1;det<=3;det++) {
    TObjArray* detArray = (TObjArray*)fArray.At(det);
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      TObjArray* vtxArray = (TObjArray*)detArray->At(ir);
      
      TH2F* hMult   = (TH2F*)vtxArray->At(vtxbin); 
      hMult->Reset();
    }
    
  }
  
  
  
  for(UShort_t det=1;det<=3;det++) {
    TObjArray* detInputArray = (TObjArray*)fInputArray->At(det);
    TObjArray* detArray = (TObjArray*)fArray.At(det);
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
      TObjArray* vtxInputArray = (TObjArray*)detInputArray->At(ir);
      TObjArray* vtxArray = (TObjArray*)detArray->At(ir);
      TH2F* hMultTotal = (TH2F*)vtxArray->At(vtxbin);
      TH2F* hMultInput = (TH2F*)vtxInputArray->At(vtxbin);
      TH2F* hHits      = (TH2F*)fOutputList->FindObject(Form("hits_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      
      hHits->Add(hMultInput);
      TH2F* hBg        = pars->GetBackgroundCorrection(det, ringChar, vtxbin);
      
      TH2F* hTmp       = (TH2F*)hMultInput->Clone("hMult_from_event");
      
      hTmp->Divide(hTmp,hBg,1,1);//,"B");
      
      hMultTotal->Add(hTmp);
      delete hTmp;
      
    }
  }
  if(fStandalone) {
    PostData(0, fOutputList); 
    PostData(1, fOutputVertexString);
  }
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBackgroundCorrection::Terminate(Option_t */*option*/) {
  
  /*
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  Int_t nVtxbins = pars->GetNvtxBins();
  
  for(UShort_t det=1;det<=3;det++) {
    TObjArray* detArray = (TObjArray*)fArray.At(det);
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      TObjArray* vtxArray = (TObjArray*)detArray->At(ir);
      for(Int_t i =0; i<nVtxbins; i++) {
	TH2F* hMultTotal = (TH2F*)vtxArray->At(i);
	if(fNevents.At(i))
	  hMultTotal->Scale(1/(Float_t)fNevents.At(i));
      }
    }
  }
  */
}
//_____________________________________________________________________
//
//
// EOF
