 
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
#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "TProfile2D.h"
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
  TH2F* hSPDMult = 0;
  TH2F* hSPDMultTrVtx = 0;
  // TH2F* hHitsNoCuts = 0;
  Int_t nVtxbins = pars->GetNvtxBins();
  for(Int_t i = 0; i< nVtxbins; i++) {
    for(Int_t det =1; det<=3;det++)
      {
	Int_t nRings = (det==1 ? 1 : 2);
	for(Int_t ring = 0;ring<nRings;ring++)
	  {
	    Char_t ringChar = (ring == 0 ? 'I' : 'O');
	    Int_t  nSec     = (ring == 0 ? 20 : 40);
	    
	    
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
	    
	  }
      }
    //HHD SPD hists
    TH2F* hBg = pars->GetBackgroundCorrection(1, 'I', i);
    hSPDMult  = new TH2F(Form("mult_SPD_vtxbin%d",i),Form("mult_SPD_vtxbin%d",i),
			 hBg->GetNbinsX(),
			 hBg->GetXaxis()->GetXmin(),
			 hBg->GetXaxis()->GetXmax(),
			 20, 0, 2*TMath::Pi());
    hSPDMult->Sumw2();
    fOutputList->Add(hSPDMult);
    hSPDMultTrVtx  = new TH2F(Form("multTrVtx_SPD_vtxbin%d",i),Form("multTrVtx_SPD_vtxbin%d",i),
			 hBg->GetNbinsX(),
			 hBg->GetXaxis()->GetXmin(),
			 hBg->GetXaxis()->GetXmax(),
			 20, 0, 2*TMath::Pi());
    hSPDMultTrVtx->Sumw2();
    fOutputList->Add(hSPDMultTrVtx);
    
  }
  
  TH2F* dNdetadphiHistogram = new TH2F("dNdetadphiHistogramTrVtx","dNdetadphiHistogramTrVtx;#eta;#Phi",pars->GetNetaBins(),-6,6,20,0,2*TMath::Pi());
  TH2F* dNdetadphiHistogramSPD = new TH2F("dNdetadphiHistogramSPDTrVtx","dNdetadphiHistogramSPDTrVtx;#eta;#Phi",pars->GetNetaBins(),-6,6,20,0,2*TMath::Pi());
  //dNdetadphiHistogram->SetErrorOption("g");
  
  fHitList->Add(dNdetadphiHistogram);
  fOutputList->Add(dNdetadphiHistogram);
  fHitList->Add(dNdetadphiHistogramSPD);
  fOutputList->Add(dNdetadphiHistogramSPD);
  
  
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
    
      TH2F* hSPDMult      = (TH2F*)fOutputList->FindObject(Form("mult_SPD_vtxbin%d",vtxbin));
      hSPDMult->Reset();
      TH2F* hSPDMultTrVtx = (TH2F*)fOutputList->FindObject(Form("multTrVtx_SPD_vtxbin%d",vtxbin));
      hSPDMultTrVtx->Reset();
    }
    
  }
  
  
  
  for(UShort_t det=1;det<=3;det++) {
    
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
      
      TH2F* hMult      = (TH2F*)fOutputList->FindObject(Form("mult_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hMultTrVtx = (TH2F*)fOutputList->FindObject(Form("multTrVtx_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hMultInput = (TH2F*)fInputList->FindObject(Form("FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      TH2F* hHits      = (TH2F*)fHitList->FindObject(Form("hits_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      
      //if(pars->GetProcessHits())
      hHits->Add(hMultInput);
      
      TH2F* hBg        = pars->GetBackgroundCorrection(det, ringChar, vtxbin);
      
      hMult->Add(hMultInput);
      hMultTrVtx->Add(hMultInput);
      hMult->Divide(hBg);//,"B");
      hMultTrVtx->Divide(hBg);//,"B");

      //sharing efficiency correction ?
      if(pars->SharingEffPresent()) {
	TH1F* hSharingEff = pars->GetSharingEfficiencyTrVtx(det,ringChar,vtxbin); 
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
      }
      
      //if(pars->GetEventSelectionEfficiency(vtxbin) > 0)
      //	hMult->Scale(1/pars->GetEventSelectionEfficiency(vtxbin));
	//else
	//hMult->Scale(0);
      hMult->Divide(pars->GetEventSelectionEfficiency(vtxbin,ringChar));
      
      
      }
  }
  
  //HHD SPD code
  
  TH2F* hSPDMult      = (TH2F*)fOutputList->FindObject(Form("mult_SPD_vtxbin%d",vtxbin));
  TH2F* hSPDMultTrVtx = (TH2F*)fOutputList->FindObject(Form("multTrVtx_SPD_vtxbin%d",vtxbin));
  
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  AliESDEvent* esd = esdH->GetEvent();
  const AliMultiplicity* spdmult = esd->GetMultiplicity();
  for(Int_t j = 0; j< spdmult->GetNumberOfTracklets();j++) {
    hSPDMult->Fill(spdmult->GetEta(j),spdmult->GetPhi(j));
    hSPDMultTrVtx->Fill(spdmult->GetEta(j),spdmult->GetPhi(j));
  }
  for(Int_t j = 0; j< spdmult->GetNumberOfSingleClusters();j++) {
    hSPDMult->Fill(-TMath::Log(TMath::Tan(spdmult->GetThetaSingle(j)/2.)),spdmult->GetPhiSingle(j));
    hSPDMultTrVtx->Fill(-TMath::Log(TMath::Tan(spdmult->GetThetaSingle(j)/2.)),spdmult->GetPhiSingle(j));
    
  }
  
  TH2F* hBgSPD        = pars->GetBackgroundCorrection(0, 'Q', vtxbin);
  if(hBgSPD) { 
  TH1F* hDead      = pars->GetSPDDeadCorrection(vtxbin);
  for(Int_t i=1; i<=hSPDMult->GetNbinsX(); i++) {
    for(Int_t j=1; j<=hSPDMult->GetNbinsY(); j++) {
      Float_t mult = hSPDMult->GetBinContent(i,j);
      Float_t correction = hBgSPD->GetBinContent(i,j);
      Float_t correctedMult = 0;
      Float_t correctedError = 0;
      if(correction > 0 && mult > 0) {
	correctedMult = mult/correction;
	if(hDead->GetBinContent(i) > 0)
	  correctedMult = correctedMult/hDead->GetBinContent(i);
	correctedError = correctedMult*TMath::Sqrt( TMath::Power(hSPDMult->GetBinError(i,j)/hSPDMult->GetBinContent(i,j),2) + 
						    TMath::Power(hBgSPD->GetBinError(i,j)/hBgSPD->GetBinContent(i,j),2));
	
      }
      
      if(correctedMult != 0) {
	hSPDMult->SetBinContent(i,j,correctedMult);
	hSPDMultTrVtx->SetBinContent(i,j,correctedMult);
	hSPDMult->SetBinError(i,j,correctedError);
	hSPDMultTrVtx->SetBinError(i,j,correctedError);
      }
    }
  }
  hSPDMult->Divide(pars->GetEventSelectionEfficiency(vtxbin,'I'));
  
  //if(pars->GetEventSelectionEfficiency(vtxbin) > 0)
  // hSPDMult->Scale(1/pars->GetEventSelectionEfficiency(vtxbin));
    //else
  //  hSPDMult->Scale(0);
  
  }
  else
    AliWarning("No SPD background map found");
  
  //std::cout<<spdmult->GetNumberOfTracklets()<<"  "<<spdmult->GetNumberOfITSClusters(0)<<"    "<< spdmult->GetNumberOfSingleClusters()<<std::endl;
  
  CreatePerEventHistogram(vtxbin);
  
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
	TH2F* hHits      = (TH2F*)fHitList->FindObject(Form("hits_FMD%d%c_vtxbin%d",det,ringChar,i));
	TH1D* hHitsproj  = hHits->ProjectionX(Form("hits_FMD%d%c_vtxbin%d_proj",det,ringChar,i),1,hHits->GetNbinsY());
	TH1D* hHitsNoCuts = (TH1D*)hHitsproj->Clone(Form("hits_NoCuts_FMD%d%c_vtxbin%d_proj",det,ringChar,i));
	if(pars->GetEventSelectionEfficiency(i) > 0)
	  hHitsNoCuts->Scale(1/pars->GetEventSelectionEfficiency(i));
	else
	  hHitsNoCuts->Scale(0);
	fHitList->Add(hHitsproj);
	fHitList->Add(hHitsNoCuts);
	
      }
    }
  }
}

//_____________________________________________________________________
void AliFMDAnalysisTaskBackgroundCorrection::CreatePerEventHistogram(Int_t vtxbin) {
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  TH2F* dNdetadphiHistogram = (TH2F*)fHitList->FindObject("dNdetadphiHistogramTrVtx");
  TH2F* dNdetadphiHistogramSPD = (TH2F*)fHitList->FindObject("dNdetadphiHistogramSPDTrVtx");
  
  dNdetadphiHistogram->Reset();
  dNdetadphiHistogramSPD->Reset();
  
  for(Int_t det = 0; det<=3; det++) {
    Int_t maxRing = (det<= 1 ? 0 : 1);
    
    
    
    for(Int_t ring = 0;ring<=maxRing;ring++) {
      
      Char_t ringChar = (ring == 0 ? 'I' : 'O');
      
      TH2F* multhistoriginal = 0;
      
      if(det == 0)
	multhistoriginal = (TH2F*)fOutputList->FindObject(Form("multTrVtx_SPD_vtxbin%d",vtxbin));
      else	 
	multhistoriginal = (TH2F*)fOutputList->FindObject(Form("multTrVtx_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
      
      TH2F* multhist = (TH2F*)multhistoriginal->Clone("tmp");
      
      
      
      if(ringChar == 'O' && det > 0)
	multhist->RebinY(2);
      
      for(Int_t i=pars->GetFirstEtaBinToInclude(vtxbin,det,ringChar); i<=pars->GetLastEtaBinToInclude(vtxbin,det,ringChar); i++) {
	for(Int_t j=1; j<=multhist->GetNbinsY(); j++) {
	  if(multhist->GetBinContent(i,j) < 0.0001) continue;
	  
	  Bool_t overlapFMD = kFALSE;
	  Bool_t overlapSPD = kFALSE;
	  
	  if(det == 1 && ringChar =='I')
	    if(i<=pars->GetLastEtaBinToInclude(vtxbin,2,'I'))
	      overlapFMD = kTRUE;
		  
	  if(det == 2 && ringChar =='O') {
	    if(i>=pars->GetFirstEtaBinToInclude(vtxbin,2,'I'))
	      overlapFMD = kTRUE;
	    if(i<=pars->GetLastEtaBinToInclude(vtxbin,0,'I'))// && TMath::Abs(multhist->GetXaxis()->GetBinCenter(i)) < 2)
	      overlapSPD = kTRUE;
	  }
	  if(det == 2 && ringChar =='I')
	    if(i<=pars->GetLastEtaBinToInclude(vtxbin,2,'O') || i>=pars->GetFirstEtaBinToInclude(vtxbin,1,'I'))
	      overlapFMD = kTRUE;
		  
	  if(det == 3 && ringChar =='I')
	    if(i>=pars->GetFirstEtaBinToInclude(vtxbin,3,'O'))
	      overlapFMD = kTRUE; 
	  
	  if(det == 3 && ringChar =='O') {
	    if(i<=pars->GetLastEtaBinToInclude(vtxbin,3,'I'))
	      overlapFMD = kTRUE;
	    if(i>=pars->GetFirstEtaBinToInclude(vtxbin,0,'I'))// && TMath::Abs(multhist->GetXaxis()->GetBinCenter(i)) < 2)
	      overlapSPD = kTRUE;
	  }
	  
	  if(det == 0) {
	    if(i<=pars->GetLastEtaBinToInclude(vtxbin,3,'O') || i>=pars->GetFirstEtaBinToInclude(vtxbin,2,'O'))
	      overlapSPD = kTRUE;
	  }
	  //std::cout<<overlapFMD<<"    "<<overlapSPD<<std::endl;
	  
	  
	  if(det > 0) {
	    if(overlapFMD) {
	      dNdetadphiHistogram->SetBinContent(i,j,dNdetadphiHistogram->GetBinContent(i,j)+0.5*multhist->GetBinContent(i,j));
	      dNdetadphiHistogram->SetBinError(i,j,TMath::Sqrt(TMath::Power(dNdetadphiHistogram->GetBinError(i,j),2)+TMath::Power(0.5*multhist->GetBinError(i,j),2)));
	    }
	    else {
	      dNdetadphiHistogram->SetBinContent(i,j,multhist->GetBinContent(i,j));
	      dNdetadphiHistogram->SetBinError(i,j,multhist->GetBinError(i,j));
	    }
	  }
	  
	  if( overlapFMD && overlapSPD) {
	    dNdetadphiHistogramSPD->SetBinContent(i,j,dNdetadphiHistogramSPD->GetBinContent(i,j)+0.33*multhist->GetBinContent(i,j));
	    dNdetadphiHistogramSPD->SetBinError(i,j,TMath::Sqrt(TMath::Power(dNdetadphiHistogramSPD->GetBinError(i,j),2)+TMath::Power(0.33*multhist->GetBinError(i,j),2)));
	  }
	  else if( overlapFMD || overlapSPD) {
	    dNdetadphiHistogramSPD->SetBinContent(i,j,dNdetadphiHistogramSPD->GetBinContent(i,j)+0.5*multhist->GetBinContent(i,j));
	    dNdetadphiHistogramSPD->SetBinError(i,j,TMath::Sqrt(TMath::Power(dNdetadphiHistogramSPD->GetBinError(i,j),2)+TMath::Power(0.5*multhist->GetBinError(i,j),2)));
	  }
	  else {
	    dNdetadphiHistogramSPD->SetBinContent(i,j,multhist->GetBinContent(i,j));
	    dNdetadphiHistogramSPD->SetBinError(i,j,multhist->GetBinError(i,j));
	  }
	  
	  
	  

	
	}
      }
      delete multhist;
    }
  }
}
//_____________________________________________________________________
//
//
// EOF
