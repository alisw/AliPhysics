 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>
#include "TH2F.h"
#include "AliFMDAnalysisTaskDndeta.h"
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

ClassImp(AliFMDAnalysisTaskDndeta)


AliFMDAnalysisTaskDndeta::AliFMDAnalysisTaskDndeta()
: fDebug(0),
  fOutputList(0),
  fInputList(0),
  fArray(0),
  fInputArray(0),
  fVertexString(0x0),
  fNevents(),
  fStandalone(kTRUE)
{
  // Default constructor
  DefineInput (0, TList::Class());
  DefineOutput(0, TList::Class());
}
//_____________________________________________________________________
AliFMDAnalysisTaskDndeta::AliFMDAnalysisTaskDndeta(const char* name, Bool_t SE):
    AliAnalysisTask(name, "Density"),
    fDebug(0),
    fOutputList(0),
    fInputList(0),
    fArray(),
    fInputArray(0),
    fVertexString(0x0),
    fNevents(),
    fStandalone(kTRUE)
{
  fStandalone = SE;
  if(fStandalone) {
    DefineInput (0, TList::Class());
    DefineInput(1, TObjString::Class());
    DefineOutput(0, TList::Class());
    
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskDndeta::CreateOutputObjects()
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  fArray.SetName("FMD");
  fArray.SetOwner();
  
  if(!fOutputList)
    fOutputList = new TList();
  fOutputList->SetName("BackgroundCorrected");
  
  
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
	    hMult  = new TH2F(Form("dNdeta_FMD%d%c_vtxbin%d",det,ringChar,i),Form("dNdeta_FMD%d%c_vtxbin%d",det,ringChar,i),
			      hBg->GetNbinsX(),
			      hBg->GetXaxis()->GetXmin(),
			      hBg->GetXaxis()->GetXmax(),
			      nSec, 0, 2*TMath::Pi());
	    hMult->Sumw2();
	    fOutputList->Add(hMult);
	    vtxArray->AddAtAndExpand(hMult,i);
	    
	  }
	} 
    }
  
  fNevents.SetBins(nVtxbins,0,nVtxbins);
  fNevents.SetName("nEvents");
  fOutputList->Add(&fNevents);
   
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskDndeta::ConnectInputData(Option_t */*option*/)
{
  if(fStandalone) {
    fInputList   = (TList*)GetInputData(0);
    fVertexString = (TObjString*)GetInputData(1);
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskDndeta::Exec(Option_t */*option*/)
{
  Int_t vtxbin   = fVertexString->GetString().Atoi();
  fNevents.Fill(vtxbin);
  
  for(UShort_t det=1;det<=3;det++) {
    //TObjArray* detInputArray = (TObjArray*)fInputArray->At(det);
    TObjArray* detArray = (TObjArray*)fArray.At(det);
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
      //TObjArray* vtxInputArray = (TObjArray*)detInputArray->At(ir);
      TObjArray* vtxArray = (TObjArray*)detArray->At(ir);
      TH2F* hMultTotal = (TH2F*)vtxArray->At(vtxbin);
      
      
      TH2F* hMultInput = (TH2F*)fInputList->FindObject(Form("mult_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      
      hMultTotal->Add(hMultInput);
      
      
    }
  }
  if(fStandalone) {
    PostData(0, fOutputList); 
  }
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskDndeta::Terminate(Option_t */*option*/) {
  
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  Int_t nVtxbins = pars->GetNvtxBins();
  
  for(UShort_t det=1;det<=3;det++) {
    TObjArray* detArray = (TObjArray*)fArray.At(det);
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      TObjArray* vtxArray = (TObjArray*)detArray->At(ir);
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
      for(Int_t i =0; i<nVtxbins; i++) {
	TH2F* hMultTotal = (TH2F*)vtxArray->At(i);
	TH1D* hMultProj   = hMultTotal->ProjectionX(Form("dNdeta_FMD%d%c_vtxbin%d_proj",det,ringChar,i),1,hMultTotal->GetNbinsY());
	fOutputList->Add(hMultProj);
      }
    }
  }
  
}
//_____________________________________________________________________
//
//
// EOF
