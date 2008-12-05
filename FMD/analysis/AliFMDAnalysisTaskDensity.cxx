 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>
#include "TAxis.h"
#include "TH2F.h"
#include "AliFMDAnalysisTaskDensity.h"
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

ClassImp(AliFMDAnalysisTaskDensity)

//_____________________________________________________________________
AliFMDAnalysisTaskDensity::AliFMDAnalysisTaskDensity()
: fDebug(0),
  fOutputList(),
  fArray(),
  fESD(0x0),
  fVertexString(),
  fVertex(0),
  fStandalone(kTRUE)
{
  // Default constructor
  DefineInput (0, AliESDFMD::Class());
  DefineInput (1, AliESDVertex::Class());
  DefineOutput(0,TList::Class());
}
//_____________________________________________________________________
AliFMDAnalysisTaskDensity::AliFMDAnalysisTaskDensity(const char* name, Bool_t SE):
    AliAnalysisTask(name, "Density"),
    fDebug(0),
    fOutputList(0),
    fArray(),
    fESD(0x0),
    fVertexString(),
    fVertex(0),
    fStandalone(kTRUE)
{
  fStandalone = SE;
  if(fStandalone) {
    DefineInput (0, AliESDFMD::Class());
    DefineInput (1, AliESDVertex::Class());
    DefineOutput(0, TList::Class());
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskDensity::CreateOutputObjects()
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  fArray.SetName("FMD");
  fArray.SetOwner();
  if(!fOutputList)
    fOutputList = new TList();
  fOutputList->SetName("density_list");
  
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
	    
	    hMult  = new TH2F(Form("FMD%d%c_vtxbin%d",det,ringChar,i),Form("FMD%d%c_vtxbin%d",det,ringChar,i),
			      hBg->GetNbinsX(),
			      hBg->GetXaxis()->GetXmin(),
			      hBg->GetXaxis()->GetXmax(),
			      nSec, 0, 2*TMath::Pi());
	    vtxArray->AddAtAndExpand(hMult,i);
	    
	  }
	} 
    }
  
  fOutputList->Add(&fArray);
  fOutputList->Add(&fVertexString);
  
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskDensity::ConnectInputData(Option_t */*option*/)
{
  if(fStandalone) {
    fESD    = (AliESDFMD*)GetInputData(0);
    fVertex = (AliESDVertex*)GetInputData(1);
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskDensity::Exec(Option_t */*option*/)
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  AliFMDGeometry* geo       = AliFMDGeometry::Instance();
  
  //AliESDFMD*   fmd = fESD->GetFMDData();
  
  Double_t vertex[3];
  fVertex->GetXYZ(vertex);
  // Z Vtx cut
  if( TMath::Abs(vertex[2]) > pars->GetVtxCutZ()) 
    return;
  Double_t delta = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
  Double_t vertexBinDouble = (vertex[2] + pars->GetVtxCutZ()) / delta;
  
  Int_t vtxbin = (Int_t)vertexBinDouble;
  
  fVertexString.SetString(Form("%d",vtxbin));
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
    TObjArray* detArray = (TObjArray*)fArray.At(det);
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      TObjArray* vtxArray = (TObjArray*)detArray->At(ir);
      
      TH2F* hMult   = (TH2F*)vtxArray->At(vtxbin);
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      for(UShort_t sec =0; sec < nsec;  sec++)  {
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  Float_t mult = fESD->Multiplicity(det,ring,sec,strip);
	  if(mult < 1 || mult == AliESDFMD::kInvalidMult) continue;
	  Float_t eta = fESD->Eta(det,ring,sec,strip);
	  Double_t x,y,z;
	  geo->Detector2XYZ(det,ring,sec,strip,x,y,z);
	  Float_t phi = TMath::ATan2(y,x);
	  if(phi<0)
	    phi = phi+2*TMath::Pi();
	  hMult->Fill(eta,phi,mult);
	  
	}
      }
    }
	
  
  }
  if(fStandalone) {
    PostData(0, fOutputList); 
  }
  
}
//_____________________________________________________________________
//
//EOF
//
