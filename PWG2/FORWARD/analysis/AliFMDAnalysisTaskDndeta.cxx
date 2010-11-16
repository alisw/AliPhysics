 
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
#include "AliGenPythiaEventHeader.h"
#include "AliHeader.h"
//#include "TDatabasePDG.h"
//#include "TParticlePDG.h"
#include "AliFMDStripIndex.h"
#include "AliESDInputHandler.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliLog.h"
ClassImp(AliFMDAnalysisTaskDndeta)


AliFMDAnalysisTaskDndeta::AliFMDAnalysisTaskDndeta()
: fDebug(0),
  fOutputList(0),
  fInputList(0),
  fVertexString(0x0),
  fNevents(),
  fNNSDevents(),
  fNMCevents(),
  fNMCNSDevents(),
  fStandalone(kTRUE),
  fLastTrackByStrip(0),
  fVtxEff(1),
  fVtxEffNSD(1)
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
    fNNSDevents(),
    fNMCevents(),
    fNMCNSDevents(),
    fStandalone(kTRUE),
    fLastTrackByStrip(0),
    fVtxEff(1),
    fVtxEffNSD(1)
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
  TH2F* hMultNSD = 0;
  TH1F* hHits = 0;
  TH2F* hMultTrVtx = 0;
  TH1F* hPrimVertexBin = 0;
  TH1F* hPrimVertexBinNSD = 0;
  
  TH2F* hBgTmp   = pars->GetBackgroundCorrection(1, 'I', 0);
  TH1F* hPrimary = new TH1F("hMultvsEta","hMultvsEta",
			    hBgTmp->GetNbinsX(),
			    hBgTmp->GetXaxis()->GetXmin(),
			    hBgTmp->GetXaxis()->GetXmax());
  hPrimary->Sumw2();
  fOutputList->Add(hPrimary);
  TH1F* hPrimaryNSD = new TH1F("hMultvsEtaNSD","hMultvsEtaNSD",
			    hBgTmp->GetNbinsX(),
			    hBgTmp->GetXaxis()->GetXmin(),
			    hBgTmp->GetXaxis()->GetXmax());
  hPrimaryNSD->Sumw2();
  fOutputList->Add(hPrimaryNSD);
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
	    hMultTrVtx  = new TH2F(Form("dNdetaTrVtx_FMD%d%c_vtxbin%d",det,ringChar,i),Form("dNdetaTrVtx_FMD%d%c_vtxbin%d",det,ringChar,i),
				   hBg->GetNbinsX(),
				   hBg->GetXaxis()->GetXmin(),
				   hBg->GetXaxis()->GetXmax(),
				   nSec, 0, 2*TMath::Pi());
	    hMultNSD  = new TH2F(Form("dNdetaNSD_FMD%d%c_vtxbin%d",det,ringChar,i),Form("dNdetaNSD_FMD%d%c_vtxbin%d",det,ringChar,i),
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
	    
	    hMultNSD->Sumw2();
	    fOutputList->Add(hMultNSD);
	    
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
    
    hPrimVertexBinNSD = new TH1F(Form("primmult_NSD_vtxbin%d",i),
				 Form("primmult_NSD_vtxbin%d",i),
				 hBgTmp->GetNbinsX(),
				 hBgTmp->GetXaxis()->GetXmin(),
				 hBgTmp->GetXaxis()->GetXmax());
    hPrimVertexBinNSD->Sumw2();
    fOutputList->Add(hPrimVertexBinNSD);
    

    //SPD part
    TH2F* hSPDMult  = new TH2F(Form("dNdeta_SPD_vtxbin%d",i),Form("dNdeta_SPD_vtxbin%d",i),
			 hBgTmp->GetNbinsX(),
			 hBgTmp->GetXaxis()->GetXmin(),
			 hBgTmp->GetXaxis()->GetXmax(),
			 20, 0, 2*TMath::Pi());
    hSPDMult->Sumw2();
    fOutputList->Add(hSPDMult);
    TH2F* hSPDMultTrVtx  = new TH2F(Form("dNdetaTrVtx_SPD_vtxbin%d",i),Form("dNdetaTrVtx_SPD_vtxbin%d",i),
			 hBgTmp->GetNbinsX(),
			 hBgTmp->GetXaxis()->GetXmin(),
			 hBgTmp->GetXaxis()->GetXmax(),
			 20, 0, 2*TMath::Pi());
    hSPDMultTrVtx->Sumw2();
    fOutputList->Add(hSPDMultTrVtx);

    TH2F* hSPDMultNSD  = new TH2F(Form("dNdetaNSD_SPD_vtxbin%d",i),Form("dNdetaNSD_SPD_vtxbin%d",i),
			 hBgTmp->GetNbinsX(),
			 hBgTmp->GetXaxis()->GetXmin(),
			 hBgTmp->GetXaxis()->GetXmax(),
			 20, 0, 2*TMath::Pi());
    hSPDMultNSD->Sumw2();
    fOutputList->Add(hSPDMultNSD);

  }
  
  //AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  //  TH2F* dNdetadphiHistogramTotal = new TH2F("dNdetadphiHistogramTotal","dNdetadphiHistogram;#eta;#Phi",pars->GetNetaBins(),-6,6,20,0,2*TMath::Pi());
  //dNdetadphiHistogramTotal->SetErrorOption("g");
  //fOutputList->Add(dNdetadphiHistogramTotal);
  
  
  
  fNevents.SetBins(nVtxbins,0,nVtxbins);
  fNevents.SetName("nEvents");
  fNNSDevents.SetBins(nVtxbins,0,nVtxbins);
  fNNSDevents.SetName("nNSDEvents");
  
  fNMCevents.SetBins(nVtxbins,0,nVtxbins);
  fNMCevents.SetName("nMCEvents");
  fNMCNSDevents.SetBins(nVtxbins,0,nVtxbins);
  fNMCNSDevents.SetName("nMCNSDEvents");
  fOutputList->Add(&fNevents);
  fOutputList->Add(&fNNSDevents);
  fOutputList->Add(&fNMCevents);
  fOutputList->Add(&fNMCNSDevents);
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
  
  //AliESDInputHandler* eventHandler = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  //AliESDEvent* esd = eventHandler->GetEvent();
  Bool_t nsd = pars->IsEventTriggered(AliFMDAnaParameters::kNSD);
  if(nsd) fNNSDevents.Fill(vtxbin);
  
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
            
      TH2F* hMultTotal      = (TH2F*)fOutputList->FindObject(Form("dNdeta_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hMultTotalTrVtx = (TH2F*)fOutputList->FindObject(Form("dNdetaTrVtx_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hMultTotalNSD   = (TH2F*)fOutputList->FindObject(Form("dNdetaNSD_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      
      
      TH2F* hMultInput      = (TH2F*)fInputList->FindObject(Form("mult_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hMultInputTrVtx = (TH2F*)fInputList->FindObject(Form("multTrVtx_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hMultInputNSD   = (TH2F*)fInputList->FindObject(Form("multNSD_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      
      hMultTotal->Add(hMultInput);
      hMultTotalTrVtx->Add(hMultInputTrVtx);
      if(nsd)
	hMultTotalNSD->Add(hMultInputNSD);
      
    }
  }
  
  //SPD code
  TH2F* hMultSPDTotal      = (TH2F*)fOutputList->FindObject(Form("dNdeta_SPD_vtxbin%d",vtxbin));
  TH2F* hMultSPDTotalTrVtx = (TH2F*)fOutputList->FindObject(Form("dNdetaTrVtx_SPD_vtxbin%d",vtxbin));
  TH2F* hMultSPDTotalNSD   = (TH2F*)fOutputList->FindObject(Form("dNdetaNSD_SPD_vtxbin%d",vtxbin));
  TH2F* hMultSPDInput      = (TH2F*)fInputList->FindObject(Form("mult_SPD_vtxbin%d",vtxbin));
  TH2F* hMultSPDInputTrVtx = (TH2F*)fInputList->FindObject(Form("multTrVtx_SPD_vtxbin%d",vtxbin));
  TH2F* hMultSPDInputNSD   = (TH2F*)fInputList->FindObject(Form("multNSD_SPD_vtxbin%d",vtxbin));
  
  hMultSPDTotal->Add(hMultSPDInput);
  hMultSPDTotalTrVtx->Add(hMultSPDInputTrVtx);
  if(nsd)
    hMultSPDTotalNSD->Add(hMultSPDInputNSD);
  
  if(pars->GetProcessPrimary())
    ProcessPrimary();
  
  //TH2F* dNdetadphiHistogram = (TH2F*)fOutputList->FindObject("dNdetadphiHistogramSPDTrVtx");
  
  //  TH2F* dNdetadphiHistogramTotal = (TH2F*)fOutputList->FindObject("dNdetadphiHistogramTotal");
  
  //  if(vtxbin == 4)
  //  dNdetadphiHistogramTotal->Add(dNdetadphiHistogram);
  
  
  
  
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
	fOutputList->Add(hMultTotal->Clone(Form("%s_orig", hMultTotal->GetName())));
	if(fVtxEff) { 
	  hMultTotal->Scale(fVtxEff);
	}
	
	TH2F* hMultTotalNSD = (TH2F*)fOutputList->FindObject(Form("dNdetaNSD_FMD%d%c_vtxbin%d",det,ringChar,i));
	fOutputList->Add(hMultTotalNSD->Clone(Form("%s_orig", hMultTotalNSD->GetName())));
	if(fVtxEffNSD) { 
	  hMultTotalNSD->Scale(fVtxEffNSD);
	}
	
	//TH2F* hMultTrVtx = (TH2F*)hMultTotal->Clone(Form("dNdeta_FMD%d%c_TrVtx_vtxbin%d",det,ringChar,i));
	TH2F* hMultTrVtx = (TH2F*)fOutputList->FindObject(Form("dNdetaTrVtx_FMD%d%c_vtxbin%d",det,ringChar,i));
	
	TH1D* hMultProj   = hMultTotal->ProjectionX(Form("dNdeta_FMD%d%c_vtxbin%d_proj",det,ringChar,i),1,hMultTotal->GetNbinsY());
	TH1D* hMultProjTrVtx   = hMultTrVtx->ProjectionX(Form("dNdeta_FMD%d%c_TrVtx_vtxbin%d_proj",det,ringChar,i),1,hMultTrVtx->GetNbinsY());
	TH1D* hMultProjNSD   = hMultTotalNSD->ProjectionX(Form("dNdetaNSD_FMD%d%c_vtxbin%d_proj",det,ringChar,i),1,hMultTotalNSD->GetNbinsY());
	
	//fOutputList->Add(hMultTrVtx);
	fOutputList->Add(hMultProj);
	fOutputList->Add(hMultProjTrVtx);
	fOutputList->Add(hMultProjNSD);
      }
    }
  }
  
  for(Int_t i =0; i<nVtxbins; i++) {
    
    TH2F* hSPDMult = (TH2F*)fOutputList->FindObject(Form("dNdeta_SPD_vtxbin%d",i));
    TH2F* hSPDMultTrVtx = (TH2F*)fOutputList->FindObject(Form("dNdetaTrVtx_SPD_vtxbin%d",i));
    TH2F* hSPDMultNSD = (TH2F*)fOutputList->FindObject(Form("dNdetaNSD_SPD_vtxbin%d",i));
    
    if(fVtxEff)
      hSPDMult->Scale(fVtxEff);
	
    if(fVtxEffNSD)
      hSPDMultNSD->Scale(fVtxEffNSD);
    
    TH1D* hMultProj   = hSPDMult->ProjectionX(Form("dNdeta_SPD_vtxbin%d_proj",i),1,hSPDMult->GetNbinsY());
    TH1D* hMultProjTrVtx   = hSPDMultTrVtx->ProjectionX(Form("dNdetaTrVtx_SPD_vtxbin%d_proj",i),1,hSPDMultTrVtx->GetNbinsY());
    TH1D* hMultProjNSD   = hSPDMultNSD->ProjectionX(Form("dNdetaNSD_SPD_vtxbin%d_proj",i),1,hSPDMultNSD->GetNbinsY());
    fOutputList->Add(hMultProj);
    fOutputList->Add(hMultProjTrVtx);
    fOutputList->Add(hMultProjNSD);
  
  }
  
  std::cout<<"FMD analysis accepted "<<fNevents.GetEntries()<<" INEL events and "<<fNNSDevents.GetEntries()<<" NSD events"<<std::endl;
}
//_____________________________________________________________________
void AliFMDAnalysisTaskDndeta::ProcessPrimary() {
  
  AliMCEventHandler* eventHandler = 
    dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()
				      ->GetMCtruthEventHandler());
  if (!eventHandler) return;

  AliMCEvent* mcEvent = eventHandler->MCEvent();
  if(!mcEvent) return;
  
  fLastTrackByStrip.Reset(-1);
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  AliMCParticle* particle = 0;
  AliStack* stack = mcEvent->Stack();
  
  TH1F* hPrimary = (TH1F*)fOutputList->FindObject("hMultvsEta");
  TH1F* hPrimaryNSD = (TH1F*)fOutputList->FindObject("hMultvsEtaNSD");
  AliHeader* header            = mcEvent->Header();
  AliGenEventHeader* genHeader = header->GenEventHeader();
  
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
  AliGenDPMjetEventHeader* dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*>(header->GenEventHeader());
  Bool_t nsd = kTRUE;
  
  if (!pythiaGenHeader && !dpmHeader) {
    std::cout<<" no pythia or dpm header!"<<std::endl;
  }
  else {
    if(pythiaGenHeader)	{
      Int_t pythiaType = pythiaGenHeader->ProcessType();
      
      if(pythiaType==92||pythiaType==93)
	nsd = kFALSE;
      
    }
    if(dpmHeader) {
      Int_t processType = dpmHeader->ProcessType();
      if(processType == 5 || processType == 6)  
	nsd = kFALSE;
    
    }
  }
  
  TArrayF vertex;
  genHeader->PrimaryVertex(vertex);
  if(TMath::Abs(vertex.At(2)) > pars->GetVtxCutZ())
    return;
  Double_t delta           = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
  Double_t vertexBinDouble = (vertex.At(2) + pars->GetVtxCutZ()) / delta;
  Int_t    vertexBin       = (Int_t)vertexBinDouble;
    
  Bool_t firstTrack = kTRUE;
  Bool_t firstTrackNSD = kTRUE;
  
  TH1F* hPrimVtxBinNSD = (TH1F*)fOutputList->FindObject(Form("primmult_NSD_vtxbin%d",vertexBin));
  TH1F* hPrimVtxBin = (TH1F*)fOutputList->FindObject(Form("primmult_vtxbin%d",vertexBin));
  
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
      hPrimVtxBin->Fill(particle->Eta());
      if(firstTrack) {
	fNMCevents.Fill(vertexBin);
	firstTrack = kFALSE;
      }
      
      if(nsd) {
	hPrimaryNSD->Fill(particle->Eta());
	hPrimVtxBinNSD->Fill(particle->Eta());
	if(firstTrackNSD) {
	  fNMCNSDevents.Fill(vertexBin);
	  firstTrackNSD = kFALSE;
	}
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
// EOF
