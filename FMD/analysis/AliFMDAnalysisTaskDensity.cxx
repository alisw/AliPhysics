 
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
#include "AliFMDParameters.h"
#include "AliFMDGeometry.h"
#include "AliFMDRing.h"

ClassImp(AliFMDAnalysisTaskDensity)

//_____________________________________________________________________
AliFMDAnalysisTaskDensity::AliFMDAnalysisTaskDensity()
: fDebug(0),
  fOutputList(),
  fArray(),
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
    fArray(),
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
  
  fFuncPos = new TF1("funcPos","pol1",0,6);
  fFuncPos->SetParameters(0.99925,0.00298301);
  fFuncNeg = new TF1("funcNeg","pol1",-6,0);
  fFuncNeg->SetParameters(0.987583,-0.0101022);

  
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
  
  fOutputList->Add(&fArray);
  fOutputList->Add(&fVertexString);
  
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
	  //Float_t eta = fESD->Eta(det,ring,sec,strip);
	  
	  if(mult == 0 || mult == AliESDFMD::kInvalidMult) continue;
	  //Particle number cut goes here...
	  //Double_t x,y,z;
	  //geo->Detector2XYZ(det,ring,sec,strip,x,y,z);
	  // Float_t phi = TMath::ATan2(y,x);
	  // if(phi<0)
	  //  phi = phi+2*TMath::Pi();
	  
	  Float_t phi = pars->GetPhiFromSector(det,ring,sec);
	  Float_t eta = pars->GetEtaFromStrip(det,ring,sec,strip,vertex[2]);
	  //std::cout<<phi<<"     "<<phicalc<<std::endl;
	  //  Float_t   r     = TMath::Sqrt(TMath::Power(x,2)+TMath::Power(y,2));
	 // Float_t   theta = TMath::ATan2(r,z-vertex[2]);
	  // Float_t   etacalc   = -1*TMath::Log(TMath::Tan(0.5*theta));
	   
	   //  std::cout<<eta<<"    "<<etacalc<<std::endl;
	   //eta = etacalc;
	     
	  Float_t m   = pars->GetMPV(det,ring,eta);
	  Float_t s   = pars->GetSigma(det,ring,eta);
	  //AliFMDParameters* recopars = AliFMDParameters::Instance();
	  
	  Float_t mult_cut = 0.15;//m-2*s;//0.15;//0.2;//m-3*s;// 0.2;//0.01;//m-2*s;//0.2;
	  if(ring == 'I')
	    mult_cut = 0.10;
	  //mult_cut = (5*recopars->GetPedestalWidth(det,ring,sec,strip))/(recopars->GetPulseGain(det,ring,sec,strip)*recopars->GetDACPerMIP());
	  
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
	  if(det == 3) 
	    correction = correction / fFuncNeg->Eval(eta);
	  else
	    correction = correction / fFuncPos->Eval(eta);
	  
	  // std::cout<<correction<<std::endl;
	  if(correction) nParticles = nParticles / correction;
	  if(nParticles > 0)
	    hMult->Fill(eta,phi,nParticles);
	  
	  //if(det == 1 && ring =='I' && nParticles >0)
	  //if(nParticles > 0)
	  //  std::cout<<det<<"    "<<ring<<"    "<<sec<<"    "<<strip<<"   "<<mult<<std::endl;
	  
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
  AliFMDRing fmdring(ring);
  fmdring.Init();
  Float_t   rad       = fmdring.GetMaxR()-fmdring.GetMinR();
  Float_t   segment   = rad / fmdring.GetNStrips();
  Float_t   radius    = fmdring.GetMinR() + segment*strip;
  
  Float_t   basearea1 = 0.5*fmdring.GetBaseStripLength(strip)*TMath::Power(radius,2);
  Float_t   basearea2 = 0.5*fmdring.GetBaseStripLength(strip)*TMath::Power((radius-segment),2);
  Float_t   basearea  = basearea1 - basearea2;
  
  Float_t   area1     = 0.5*fmdring.GetStripLength(strip)*TMath::Power(radius,2);
  Float_t   area2     = 0.5*fmdring.GetStripLength(strip)*TMath::Power((radius-segment),2);
  Float_t   area      = area1 - area2;
  
  Float_t correction = area/basearea;
  
  return correction;
}
//_____________________________________________________________________
Float_t AliFMDAnalysisTaskDensity::GetPhiFromSector(UShort_t det, Char_t ring, UShort_t sec)
{
  Int_t nsec = (ring == 'I' ? 20 : 40);
  Float_t basephi = 0;
  if(det == 1) 
    basephi = 1.72787594; 
  if(det == 2 && ring == 'I')
    basephi = 0.15707963;
  if(det == 2 && ring == 'O')
    basephi = 0.078539818;
  if(det == 3 && ring == 'I')
    basephi = 2.984513044;
  if(det == 3 && ring == 'O')
    basephi = 3.06305289;
  
  Float_t step = 2*TMath::Pi() / nsec;
  Float_t phi = 0;
  if(det == 3)
    phi = basephi - sec*step;
  else
    phi = basephi + sec*step;
  
  if(phi < 0) 
    phi = phi +2*TMath::Pi();
  if(phi > 2*TMath::Pi() )
    phi = phi - 2*TMath::Pi();
  
  return phi;
}


//
//EOF
//
