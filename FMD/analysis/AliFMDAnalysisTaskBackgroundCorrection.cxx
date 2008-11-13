 
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
#include "AliESDVertex.h"
#include "TMath.h"
#include "AliFMDAnaParameters.h"
#include "AliFMDGeometry.h"

ClassImp(AliFMDAnalysisTaskBackgroundCorrection)


AliFMDAnalysisTaskBackgroundCorrection::AliFMDAnalysisTaskBackgroundCorrection()
: fDebug(0),
  fOutputList(0),
  fArray(0),
  fInputArray(0),
  fVertexString(0x0),
  fNevents()
{
  // Default constructor
  DefineInput (0, TList::Class());
  DefineOutput(0, TList::Class());
}
//_____________________________________________________________________
AliFMDAnalysisTaskBackgroundCorrection::AliFMDAnalysisTaskBackgroundCorrection(const char* name):
    AliAnalysisTask(name, "Density"),
    fDebug(0),
    fOutputList(),
    fArray(),
    fInputArray(0),
    fVertexString(0x0),
    fNevents()
{
  DefineInput (0, TList::Class());
  DefineOutput(0, TList::Class());
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBackgroundCorrection::CreateOutputObjects()
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  fArray.SetName("FMD");
  fArray.SetOwner();
  
  TH2F* hMult = 0;
  
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
	    fOutputList.Add(hMult);
	    vtxArray->AddAtAndExpand(hMult,i);
	    
	  }
	} 
    }
  fNevents.Set(nVtxbins);
  //  fOutputList.Add(&fArray);
   
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBackgroundCorrection::ConnectInputData(Option_t */*option*/)
{

  TList* list = (TList*)GetInputData(0);
  fInputArray = (TObjArray*)list->At(0);
  fVertexString = (TObjString*)list->At(1);
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBackgroundCorrection::Exec(Option_t */*option*/)
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  Int_t vtxbin   = fVertexString->GetString().Atoi();
  
  fNevents.operator[](vtxbin)++;
  
  for(UShort_t det=1;det<=3;det++) {
    TObjArray* detInputArray = (TObjArray*)fInputArray->At(det);
    TObjArray* detArray = (TObjArray*)fArray.At(det);
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
      TObjArray* vtxInputArray = (TObjArray*)detInputArray->At(ir);
      TObjArray* vtxArray = (TObjArray*)detArray->At(ir);
      TH2F* hMultTotal = (TH2F*)vtxArray->At(vtxbin);
      TH2F* hMult      = (TH2F*)vtxInputArray->At(vtxbin);
      TH2F* hBg        = pars->GetBackgroundCorrection(det, ringChar, vtxbin);
      TH2F* hTmp       = (TH2F*)hMult->Clone("hMult_from_event");
      hTmp->Divide(hTmp,hBg,1,1,"B");
      
      hMultTotal->Add(hTmp);
      delete hTmp;
      
    }
  }
  
  PostData(0, &fOutputList); 
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBackgroundCorrection::Terminate(Option_t */*option*/) {
  
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
  
}
//_____________________________________________________________________
//
//
// EOF
