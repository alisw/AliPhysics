 
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
  fVertexString(0x0),
  fNevents(),
  fNMCevents(),
  fStandalone(kTRUE),
  fLastTrackByStrip(0)
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
    fVertexString(0x0),
    fNevents(),
    fNMCevents(),
    fStandalone(kTRUE),
    fLastTrackByStrip(0)
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
  
  if(!fOutputList)
    fOutputList = new TList();
  fOutputList->SetName("BackgroundCorrected");
  
  
  TH2F* hMult = 0;
  TH1F* hHits = 0;
  TH2F* hMultTrVtx = 0;
  TH1F* hPrimVertexBin = 0;
  
  
  TH2F* hBgTmp   = pars->GetBackgroundCorrection(1, 'I', 0);
  TH1F* hPrimary = new TH1F("hMultvsEta","hMultvsEta",
			    hBgTmp->GetNbinsX(),
			    hBgTmp->GetXaxis()->GetXmin(),
			    hBgTmp->GetXaxis()->GetXmax());
  hPrimary->Sumw2();
  fOutputList->Add(hPrimary);
  Int_t nVtxbins = pars->GetNvtxBins();
  TH2F* hBg = 0;
  for(Int_t i = 0; i< nVtxbins; i++) {
     
    for(Int_t det =1; det<=3;det++)
      {
	Int_t nRings = (det==1 ? 1 : 2);
	for(Int_t ring = 0;ring<nRings;ring++)
	  {
	    Char_t ringChar = (ring == 0 ? 'I' : 'O');
	    Int_t  nSec     = (ring == 0 ? 20 : 40);
	    
	    
	    
	    hBg = pars->GetBackgroundCorrection(det, ringChar, i);
	    hMult  = new TH2F(Form("dNdeta_FMD%d%c_vtxbin%d",det,ringChar,i),Form("dNdeta_FMD%d%c_vtxbin%d",det,ringChar,i),
			      hBg->GetNbinsX(),
			      hBg->GetXaxis()->GetXmin(),
			      hBg->GetXaxis()->GetXmax(),
			      nSec, 0, 2*TMath::Pi());
	    hMultTrVtx  = new TH2F(Form("dNdetaTrVtx_FMD%d%c_vtxbin%d",det,ringChar,i),Form("dNdeta_FMD%d%c_vtxbin%d",det,ringChar,i),
				   hBg->GetNbinsX(),
				   hBg->GetXaxis()->GetXmin(),
				   hBg->GetXaxis()->GetXmax(),
				   nSec, 0, 2*TMath::Pi());
	    hHits  = new TH1F(Form("hMCHits_FMD%d%c_vtxbin%d",det,ringChar,i),Form("hMCHits_FMD%d%c_vtxbin%d",det,ringChar,i),
			      hBg->GetNbinsX(),
			      hBg->GetXaxis()->GetXmin(),
			      hBg->GetXaxis()->GetXmax());
	    hHits->Sumw2();
	    fOutputList->Add(hHits);
	    
	    hMult->Sumw2();
	    fOutputList->Add(hMult);

	    hMultTrVtx->Sumw2();
	    fOutputList->Add(hMultTrVtx);
	    
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
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskDndeta::Exec(Option_t */*option*/)
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  fVertexString = (TObjString*)fInputList->At(0);
  
  Int_t vtxbin   = fVertexString->GetString().Atoi();
  
  fNevents.Fill(vtxbin);
  
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
            
      TH2F* hMultTotal = (TH2F*)fOutputList->FindObject(Form("dNdeta_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hMultTotalTrVtx = (TH2F*)fOutputList->FindObject(Form("dNdetaTrVtx_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      
      TH2F* hMultInput = (TH2F*)fInputList->FindObject(Form("mult_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hMultInputTrVtx = (TH2F*)fInputList->FindObject(Form("multTrVtx_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      
      hMultTotal->Add(hMultInput);
      hMultTotalTrVtx->Add(hMultInputTrVtx);
            
    }
  }
  
  if(pars->GetProcessPrimary())
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
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
      for(Int_t i =0; i<nVtxbins; i++) {
	
	TH2F* hMultTotal = (TH2F*)fOutputList->FindObject(Form("dNdeta_FMD%d%c_vtxbin%d",det,ringChar,i));
	//TH2F* hMultTrVtx = (TH2F*)hMultTotal->Clone(Form("dNdeta_FMD%d%c_TrVtx_vtxbin%d",det,ringChar,i));
	TH2F* hMultTrVtx = (TH2F*)fOutputList->FindObject(Form("dNdetaTrVtx_FMD%d%c_vtxbin%d",det,ringChar,i));
	/*
	TH1F* hSharingEff = pars->GetSharingEfficiency(det,ringChar,i);
	TH1F* hSharingEffTrVtx = pars->GetSharingEfficiencyTrVtx(det,ringChar,i);	
	
	
	for(Int_t nx=1; nx<hMultTotal->GetNbinsX(); nx++) {
	  Float_t correction = hSharingEff->GetBinContent(nx);
	  Float_t correctionTrVtx = hSharingEffTrVtx->GetBinContent(nx);
	  for(Int_t ny=1; ny<hMultTotal->GetNbinsY(); ny++) {
	    
	    if(correction != 0){
	      hMultTotal->SetBinContent(nx,ny,hMultTotal->GetBinContent(nx,ny)/correction);
	      Float_t error = TMath::Sqrt(TMath::Power(hMultTotal->GetBinError(nx,ny),2) + TMath::Power(hMultTotal->GetBinContent(nx,ny)*hSharingEff->GetBinError(nx),2)) / correction;
	      hMultTotal->SetBinError(nx,ny,error);
	    }
	    if(correctionTrVtx != 0){
	      hMultTrVtx->SetBinContent(nx,ny,hMultTrVtx->GetBinContent(nx,ny)/correctionTrVtx);
	      Float_t error = TMath::Sqrt(TMath::Power(hMultTrVtx->GetBinError(nx,ny),2) + TMath::Power(hMultTrVtx->GetBinContent(nx,ny)*hSharingEffTrVtx->GetBinError(nx),2)) / correctionTrVtx;
	      hMultTrVtx->SetBinError(nx,ny,error);
	    }
	  }
	  
	}
	
	//hMultTotal->Divide(hSharingEff);
	
	hMultTotal->Scale(1/pars->GetEventSelectionEfficiency(i));
	*/
	//hMultTrVtx->Divide(hSharingEffTrVtx);
	
	TH1D* hMultProj   = hMultTotal->ProjectionX(Form("dNdeta_FMD%d%c_vtxbin%d_proj",det,ringChar,i),1,hMultTotal->GetNbinsY());
	TH1D* hMultProjTrVtx   = hMultTrVtx->ProjectionX(Form("dNdeta_FMD%d%c_TrVtx_vtxbin%d_proj",det,ringChar,i),1,hMultTotal->GetNbinsY());
	fOutputList->Add(hMultTrVtx);
	fOutputList->Add(hMultProj);
	fOutputList->Add(hMultProjTrVtx);
      }
    }
  }
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskDndeta::ProcessPrimary() {
  
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  AliMCEvent* mcEvent = eventHandler->MCEvent();
  if(!mcEvent)
    return;
  
  fLastTrackByStrip.Reset(-1);
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  AliMCParticle* particle = 0;
  AliStack* stack = mcEvent->Stack();
  
  TH1F* hPrimary = (TH1F*)fOutputList->FindObject("hMultvsEta");
  AliHeader* header            = mcEvent->Header();
  AliGenEventHeader* genHeader = header->GenEventHeader();
  
  TArrayF vertex;
  genHeader->PrimaryVertex(vertex);
  if(TMath::Abs(vertex.At(2)) > pars->GetVtxCutZ())
    return;
  Double_t delta           = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
  Double_t vertexBinDouble = (vertex.At(2) + pars->GetVtxCutZ()) / delta;
  Int_t    vertexBin       = (Int_t)vertexBinDouble;
    
  Bool_t firstTrack = kTRUE;
  
  // we loop over the primaries only unless we need the hits (diagnostics running slowly)
  Int_t nTracks = stack->GetNprimary();
  if(pars->GetProcessHits())
    nTracks = stack->GetNtrack();
  
  for(Int_t i = 0 ;i<nTracks;i++) {
    particle = (AliMCParticle*) mcEvent->GetTrack(i);
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
    if(pars->GetProcessHits()) {
           
      for(Int_t j=0; j<particle->GetNumberOfTrackReferences();j++) {
	
	AliTrackReference* ref = particle->GetTrackReference(j);
	UShort_t det,sec,strip;
	Char_t   ring;
	if(ref->DetectorId() != AliTrackReference::kFMD)
	  continue;
	AliFMDStripIndex::Unpack(ref->UserId(),det,ring,sec,strip);
	Float_t thisStripTrack = fLastTrackByStrip(det,ring,sec,strip);
	if(particle->Charge() != 0 && i != thisStripTrack ) {
	  //Double_t x,y,z;
	  
	  Float_t   eta   = pars->GetEtaFromStrip(det,ring,sec,strip,vertex.At(2));//-1*TMath::Log(TMath::Tan(0.5*theta));
	  TH1F* hHits = (TH1F*)fOutputList->FindObject(Form("hMCHits_FMD%d%c_vtxbin%d",det,ring,vertexBin));
	  
	
	  hHits->Fill(eta);
	  
	  Float_t nstrips = (ring =='O' ? 256 : 512);
	  
	  fLastTrackByStrip(det,ring,sec,strip) = (Float_t)i;
	
	  if(strip >0)
	    fLastTrackByStrip(det,ring,sec,strip-1) = (Float_t)i;
	  if(strip < (nstrips - 1))
	    fLastTrackByStrip(det,ring,sec,strip+1) = (Float_t)i;
	  
	}
      
	
      }
      
      
    }

      
      
  }
}
//_____________________________________________________________________
//
//
// EOF
