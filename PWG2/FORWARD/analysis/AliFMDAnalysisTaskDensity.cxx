 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>
#include "TAxis.h"
#include "TH2F.h"
#include "TF1.h"
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
//#include "AliFMDParameters.h"
//#include "AliFMDGeometry.h"
//#include "AliFMDRing.h"

ClassImp(AliFMDAnalysisTaskDensity)

//_____________________________________________________________________
AliFMDAnalysisTaskDensity::AliFMDAnalysisTaskDensity()
: fDebug(0),
  fOutputList(),
  fESD(0x0),
  fVertexString(),
  fVertex(0),
  fStandalone(kTRUE),
  fStatus(kTRUE)
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
    fESD(0x0),
    fVertexString(),
    fVertex(0),
    fStandalone(kTRUE),
    fStatus(kTRUE)
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
  
  if(!fOutputList)
    fOutputList = new TList();
  fOutputList->SetName("density_list");
  
  fOutputList->Add(&fVertexString);
  
  
  
  TH2F* hMult = 0;
  
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
	    
	    hMult  = new TH2F(Form("FMD%d%c_vtxbin%d",det,ringChar,i),Form("FMD%d%c_vtxbin%d",det,ringChar,i),
			      hBg->GetNbinsX(),
			      hBg->GetXaxis()->GetXmin(),
			      hBg->GetXaxis()->GetXmax(),
			      nSec, 0, 2*TMath::Pi());
	    
	    fOutputList->Add(hMult);
	  }
	} 
    }
  
  
  
  
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
  // AliFMDGeometry* geo       = AliFMDGeometry::Instance();
  
  //AliESDFMD*   fmd = fESD->GetFMDData();
  
  Double_t vertex[3];
  fVertex->GetXYZ(vertex);
  // Z Vtx cut
  if( TMath::Abs(vertex[2]) > pars->GetVtxCutZ()) {
    fStatus = kFALSE;
    return;
  }
  else
    fStatus = kTRUE;
  
  Double_t delta = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
  Double_t vertexBinDouble = (vertex[2] + pars->GetVtxCutZ()) / delta;
  
  Int_t vtxbin = (Int_t)vertexBinDouble;
  
  fVertexString.SetString(Form("%d",vtxbin));
  //Reset everything
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      TH2F* hMult   = (TH2F*)fOutputList->FindObject(Form("FMD%d%c_vtxbin%d",det,ring,vtxbin));
      hMult->Reset();
    }
    
  }
  
  
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      TH2F* hMult   = (TH2F*)fOutputList->FindObject(Form("FMD%d%c_vtxbin%d",det,ring,vtxbin));
     
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      
      for(UShort_t sec =0; sec < nsec;  sec++)  {
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  Float_t mult = fESD->Multiplicity(det,ring,sec,strip);
		  
	  if(mult == 0 || mult == AliESDFMD::kInvalidMult) continue;
			  
	  Float_t phi = pars->GetPhiFromSector(det,ring,sec);
	  Float_t eta = pars->GetEtaFromStrip(det,ring,sec,strip,vertex[2]);
	  
	  Float_t mult_cut = 0.15;//m-2*s;//0.15;//0.2;//m-3*s;// 0.2;//0.01;//m-2*s;//0.2;
	    if(ring == 'I')
	      mult_cut = 0.10;
	  //Float_t mult_cut = pars->GetMPV(det,ring,eta) - 5*pars->GetSigma(det,ring,eta);
	  Float_t nParticles = 0;
	  if(fESD->GetUniqueID() == kTRUE) {
	    //proton + proton
	    
	    if(mult > mult_cut) {
	      nParticles = 1; 
	    
	    }
	  }
	  else {
	    
	    //Pb+Pb
	    Float_t mpv   = pars->GetMPV(det,ring,eta);
	    Float_t sigma = pars->GetSigma(det,ring,eta);
	    Float_t alpha = pars->Get2MIPWeight(det,ring,eta);
	    Float_t beta  = pars->Get3MIPWeight(det,ring,eta);
	    
	    Float_t sumCor = TMath::Landau(mult,mpv,sigma,kTRUE)+
	      alpha*TMath::Landau(mult,2*mpv+2*sigma*TMath::Log(2),2*sigma,kTRUE)+
	      beta*TMath::Landau(mult,3*mpv+3*sigma*TMath::Log(3),3*sigma,kTRUE);
	    Float_t weight = TMath::Landau(mult,mpv,sigma,kTRUE)+
	      2*alpha*TMath::Landau(mult,2*mpv+2*sigma*TMath::Log(2),2*sigma,kTRUE)+
	      3*beta*TMath::Landau(mult,3*mpv+3*sigma*TMath::Log(3),3*sigma,kTRUE);
	    
	    
	    if(mult > mult_cut) {
	      if(sumCor) nParticles = weight / sumCor;
	      else nParticles = 1;
	      
	    }
	    //std::cout<<sumCor<<"    "<<weight<<"    "<<"    "<<mult<<"  "<<nParticles<<std::endl;
	    
	  }
	  
	  
	  
	  
	  
	  Float_t correction = GetAcceptanceCorrection(ring,strip);
	  
	  //std::cout<<"before "<<correction<<std::endl;
	  if(fESD->GetUniqueID() == kTRUE) {
	    TH1F* hDoubleHitCorrection = pars->GetDoubleHitCorrection(det,ring);
	    
	    if(hDoubleHitCorrection->GetBinContent(hDoubleHitCorrection->FindBin(eta)) != 0)
	      correction = correction*hDoubleHitCorrection->GetBinContent(hDoubleHitCorrection->FindBin(eta));
	    
	  }
	  
	  if(correction) nParticles = nParticles / correction;
	  if(nParticles > 0)
	    hMult->Fill(eta,phi,nParticles);
	  
	  
	}
      }
      
    }
    
	
  
  }
  

  if(fStandalone) {
    PostData(0, fOutputList); 
  }
  
}
//_____________________________________________________________________
Float_t AliFMDAnalysisTaskDensity::GetAcceptanceCorrection(Char_t ring, UShort_t strip)
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  //AliFMDRing fmdring(ring);
  //fmdring.Init();
  Float_t   rad       = pars->GetMaxR(ring)-pars->GetMinR(ring);
  Float_t   nstrips   = (ring == 'I' ? 512 : 256);
  Float_t   segment   = rad / nstrips;
  Float_t   radius    = pars->GetMinR(ring) + segment*strip;
  
  Float_t   basearea1 = 0.5*pars->GetBaseStripLength(ring,strip)*TMath::Power(radius,2);
  Float_t   basearea2 = 0.5*pars->GetBaseStripLength(ring,strip)*TMath::Power((radius-segment),2);
  Float_t   basearea  = basearea1 - basearea2;
  
  Float_t   area1     = 0.5*pars->GetStripLength(ring,strip)*TMath::Power(radius,2);
  Float_t   area2     = 0.5*pars->GetStripLength(ring,strip)*TMath::Power((radius-segment),2);
  Float_t   area      = area1 - area2;
  
  Float_t correction = area/basearea;
  
  return correction;
}
//_____________________________________________________________________
//
//EOF
//
