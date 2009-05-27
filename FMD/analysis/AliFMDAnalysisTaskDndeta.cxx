 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>
#include "TH1F.h"
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
//#include "AliFMDGeometry.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
//#include "TDatabasePDG.h"
//#include "TParticlePDG.h"
#include "AliFMDStripIndex.h"
ClassImp(AliFMDAnalysisTaskDndeta)


AliFMDAnalysisTaskDndeta::AliFMDAnalysisTaskDndeta()
: fDebug(0),
  fOutputList(0),
  fInputList(0),
  fArray(0),
  fInputArray(0),
  fVertexString(0x0),
  fNevents(),
  fNMCevents(),
  fStandalone(kTRUE),
  fMCevent(0),
  fLastTrackByStrip()
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
    fNMCevents(),
    fStandalone(kTRUE),
    fMCevent(0),
    fLastTrackByStrip()
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
  TH1F* hHits = 0;
  TH1F* hPrimVertexBin = 0;
  
  
  TH2F* hBgTmp   = pars->GetBackgroundCorrection(1, 'I', 0);
  TH1F* hPrimary = new TH1F("hMultvsEta","hMultvsEta",
			    hBgTmp->GetNbinsX(),
			    hBgTmp->GetXaxis()->GetXmin(),
			    hBgTmp->GetXaxis()->GetXmax());
  hPrimary->Sumw2();
  fOutputList->Add(hPrimary);
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
	    
	    hHits  = new TH1F(Form("hHits_FMD%d%c_vtxbin%d",det,ringChar,i),Form("hHits_FMD%d%c_vtxbin%d",det,ringChar,i),
			      hBg->GetNbinsX(),
			      hBg->GetXaxis()->GetXmin(),
			      hBg->GetXaxis()->GetXmax());
			    
	    
	    
	    hMult->Sumw2();
	    hHits->Sumw2();
	    fOutputList->Add(hMult);
	    fOutputList->Add(hHits);
	    vtxArray->AddAtAndExpand(hMult,i);
	    
	  }
	} 
    }
  
  for(Int_t i = 0; i< nVtxbins; i++) {
   
    hPrimVertexBin = new TH1F(Form("primmult_vtxbin%d",i),
			      Form("primmult_vtxbin%d",i),
			      hBgTmp->GetNbinsX(),
			      hBgTmp->GetXaxis()->GetXmin(),
			      hBgTmp->GetXaxis()->GetXmax());
    hPrimVertexBin->Sumw2();
    fOutputList->Add(hPrimVertexBin);
    
  }
  
  fNevents.SetBins(nVtxbins,0,nVtxbins);
  fNevents.SetName("nEvents");
  fNMCevents.SetBins(nVtxbins,0,nVtxbins);
  fNMCevents.SetName("nMCEvents");
  fOutputList->Add(&fNevents);
  fOutputList->Add(&fNMCevents);
  
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
  
  if(fMCevent)
    ProcessPrimary();
  
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
void AliFMDAnalysisTaskDndeta::ProcessPrimary() {
  
  fLastTrackByStrip.Reset(-1);
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  AliMCParticle* particle = 0;
  AliStack* stack = fMCevent->Stack();
  
  TH1F* hPrimary = (TH1F*)fOutputList->FindObject("hMultvsEta");
  AliHeader* header            = fMCevent->Header();
  AliGenEventHeader* genHeader = header->GenEventHeader();
  
  TArrayF vertex;
  genHeader->PrimaryVertex(vertex);
  if(TMath::Abs(vertex.At(2)) > pars->GetVtxCutZ())
    return;
  Double_t delta           = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
  Double_t vertexBinDouble = (vertex.At(2) + pars->GetVtxCutZ()) / delta;
  Int_t    vertexBin       = (Int_t)vertexBinDouble;
  
  Bool_t firstTrack = kTRUE;
  Int_t nTracks = fMCevent->GetNumberOfTracks();
  for(Int_t i = 0 ;i<nTracks;i++) {
    particle = fMCevent->GetTrack(i);
    if(!particle)
      continue;
    
    if(stack->IsPhysicalPrimary(i) && particle->Charge() != 0) {
      hPrimary->Fill(particle->Eta());
      

      TH1F* hPrimVtxBin = (TH1F*)fOutputList->FindObject(Form("primmult_vtxbin%d",vertexBin));
      hPrimVtxBin->Fill(particle->Eta());
      if(firstTrack) {
	fNMCevents.Fill(vertexBin);
	firstTrack = kFALSE;
      }
    }
     
    for(Int_t j=0; j<particle->GetNumberOfTrackReferences();j++) {
      
      AliTrackReference* ref = particle->GetTrackReference(j);
      UShort_t det,sec,strip;
      Char_t   ring;
      if(ref->DetectorId() != AliTrackReference::kFMD)
	continue;
      AliFMDStripIndex::Unpack(ref->UserId(),det,ring,sec,strip);
      Float_t thisStripTrack = fLastTrackByStrip.operator()(det,ring,sec,strip);
      if(particle->Charge() != 0 && i != thisStripTrack ) {
	//Double_t x,y,z;
	/*AliFMDGeometry* fmdgeo = AliFMDGeometry::Instance();
	fmdgeo->Detector2XYZ(det,ring,sec,strip,x,y,z);
	
	Float_t   phi   = TMath::ATan2(y,x);
	if(phi<0) phi   = phi+2*TMath::Pi();
	Float_t   r     = TMath::Sqrt(TMath::Power(x,2)+TMath::Power(y,2));
	Float_t   theta = TMath::ATan2(r,z-vertex.At(2));*/
	Float_t   eta   = pars->GetEtaFromStrip(det,ring,sec,strip,vertex.At(2));//-1*TMath::Log(TMath::Tan(0.5*theta));
	TH1F* hHits = (TH1F*)fOutputList->FindObject(Form("hHits_FMD%d%c_vtxbin%d",det,ring,vertexBin));
	hHits->Fill(eta);
	Float_t nstrips = (ring =='O' ? 256 : 512);
	
	//if(det == 1 && ring == 'I')
	//	std::cout<<"hit in "<<det<<"   "<<ring<<"   "<<sec<<"   "<<strip<<"   "<<std::endl;
	fLastTrackByStrip.operator()(det,ring,sec,strip) = (Float_t)i;
	
	if(strip >0)
	  fLastTrackByStrip.operator()(det,ring,sec,strip-1) = (Float_t)i;
	if(strip < (nstrips - 1))
	  fLastTrackByStrip.operator()(det,ring,sec,strip+1) = (Float_t)i;
	
	
      }
    }
    
    
  }
  
}
//_____________________________________________________________________
//
//
// EOF
