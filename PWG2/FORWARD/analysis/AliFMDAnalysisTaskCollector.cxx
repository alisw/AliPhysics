 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TMath.h>
#include <iostream>

#include "AliFMDAnalysisTaskCollector.h"
#include "AliAnalysisManager.h"
#include "AliESDFMD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliESDVertex.h"
#include "AliFMDAnaParameters.h"
//#include "AliFMDGeometry.h"


ClassImp(AliFMDAnalysisTaskCollector)

//____________________________________________________________________
Double_t  AliFMDAnalysisTaskCollector::TripleLandau(Double_t *x, Double_t *par) {
  
  Double_t energy        = x[0];
  Double_t constant = par[0];
  Double_t mpv      = par[1];
  Double_t sigma    = par[2];
  Double_t alpha    = par[3];
  Double_t beta     = par[4];
 
  Double_t f = constant*(TMath::Landau(energy,mpv,sigma,kTRUE)+
			 alpha*TMath::Landau(energy,2*mpv+2*sigma*TMath::Log(2),2*sigma,kTRUE)+
			 beta*TMath::Landau(energy,3*mpv+3*sigma*TMath::Log(3),3*sigma,kTRUE) );
  
  return f;
}
//____________________________________________________________________

AliFMDAnalysisTaskCollector::AliFMDAnalysisTaskCollector()
: fDebug(0),
  fOutputList(0),
  fArray(0),
  fZvtxDist(0)
{
  // Default constructor
  
 
}
//____________________________________________________________________
AliFMDAnalysisTaskCollector::AliFMDAnalysisTaskCollector(const char* name):
    AliAnalysisTaskSE(name),
    fDebug(0),
    fOutputList(0),
    fArray(0),
    fZvtxDist(0)
{
  // Default constructor
  
  DefineOutput(1, TList::Class());
}
//____________________________________________________________________
void AliFMDAnalysisTaskCollector::UserCreateOutputObjects()
{
  // Create the output container
  printf("AnalysisTaskFMD::CreateOutPutData() \n");
  
  fOutputList = new TList();//(TList*)GetOutputData(0);
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  fArray     = new TObjArray();
  fArray->SetName("FMD");
  fArray->SetOwner();
  TH1F* hEdist = 0;
 
  for(Int_t nEta = 0; nEta <= pars->GetNetaBins()+1; nEta++) {
    TObjArray* etaArray = new TObjArray();
    fArray->AddAtAndExpand(etaArray,nEta);
    for(Int_t det =1; det<=3;det++)
      {
	TObjArray* detArray = new TObjArray();
	detArray->SetName(Form("FMD%d",det));
	etaArray->AddAtAndExpand(detArray,det);
	Int_t nRings = (det==1 ? 1 : 2);
	for(Int_t ring = 0;ring<nRings;ring++)
	  {
	    Char_t ringChar = (ring == 0 ? 'I' : 'O');
	    hEdist = new TH1F(Form("FMD%d%c_etabin%d",det,ringChar,nEta),Form("FMD%d%c_etabin%d",det,ringChar,nEta),200,0,6);
	    hEdist->SetXTitle("#Delta E / E_{MIP}");
	    fOutputList->Add(hEdist);
	    detArray->AddAtAndExpand(hEdist,ring);
	  } 
      }
    
  }
  
  fZvtxDist  = new TH1F("ZvtxDist","Vertex distribution",100,-30,30);
  fZvtxDist->SetXTitle("z vertex");
  fOutputList->Add(fZvtxDist);
}

//____________________________________________________________________
void AliFMDAnalysisTaskCollector::UserExec(Option_t */*option*/)
{
  
  AliESDEvent* esd = (AliESDEvent*)InputEvent();
  AliESD* old = esd->GetAliESDOld();
  if (old) {
    esd->CopyFromOldESD();
  }
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  Bool_t trigger = pars->IsEventTriggered(esd);
  if(!trigger)
    return;
  Double_t vertex[3];
  
  pars->GetVertex(esd,vertex);
  if(vertex[0] == 0 && vertex[1] == 0 && vertex[2] == 0)
    return;
  
  fZvtxDist->Fill(vertex[2]);
  
  if(TMath::Abs(vertex[2]) > pars->GetVtxCutZ())
    return;
  
  AliESDFMD* fmd = esd->GetFMDData();
  if (!fmd) return;
  
  for(UShort_t det=1;det<=3;det++) {
      
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
  
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      TH2F* hBg = pars->GetBackgroundCorrection(det,ring,0);
      
      for(UShort_t sec =0; sec < nsec;  sec++)  {
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  
	  
	  Float_t mult = fmd->Multiplicity(det,ring,sec,strip);
	  if(mult == AliESDFMD::kInvalidMult || mult == 0) continue;
	  
	  Float_t eta = pars->GetEtaFromStrip(det,ring,sec,strip,vertex[2]);
	  
	  Int_t nEta = hBg->GetXaxis()->FindBin(eta);
	  
	  TObjArray* etaArray = (TObjArray*)fArray->At(nEta);
	  TObjArray* detArray = (TObjArray*)etaArray->At(det);
	  TH1F* Edist = (TH1F*)detArray->At(ir);
	  
	  Edist->Fill(mult);
	  
	}
      }
    }
  }
  
  PostData(1, fOutputList); 
  
}

//____________________________________________________________________
void AliFMDAnalysisTaskCollector::ReadFromFile(const Char_t* filename, Bool_t store, Int_t speciesOption)  {

  //speciesOption:
  //0: p+p Landau fit
  //1: Pb+Pb triple landau convolution fit
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  pars->Init(kTRUE,AliFMDAnaParameters::kBackgroundCorrection);
  
  TFile fin(filename,"UPDATE");
  
  TList* list = (TList*)fin.Get("energyDist");
  
  AliFMDAnaCalibEnergyDistribution* EnergyDist = new AliFMDAnaCalibEnergyDistribution();
  
  EnergyDist->SetNetaBins(pars->GetNetaBins());
  EnergyDist->SetEtaLimits(pars->GetEtaMin(),pars->GetEtaMax());
    
  for(Int_t nEta = 1; nEta <= pars->GetNetaBins(); nEta++) {
  
    for(Int_t det = 1; det<=3; det++) {
      Int_t nRings  =  (det == 1 ? 1 : 2);
      for(Int_t ring = 0;ring<nRings; ring++) {
	Char_t ringChar = (ring == 0 ? 'I' : 'O');
	
	TH1F* hEdist = (TH1F*)list->FindObject(Form("FMD%d%c_etabin%d",det,ringChar,nEta));
	TF1* fitFunc = 0 ;
	
	if(hEdist->GetEntries() != 0) {
	  
	  hEdist->GetXaxis()->SetRangeUser(0.2,hEdist->GetXaxis()->GetXmax());
	  
	  if(speciesOption == 0)
	    fitFunc =  new TF1("FMDfitFunc","landau",hEdist->GetBinCenter(hEdist->GetMaximumBin())-0.2,3);
	  if(speciesOption == 1) {
	    fitFunc = new TF1("FMDfitFunc",TripleLandau,hEdist->GetBinCenter(hEdist->GetMaximumBin())-0.2,5,5);
	    fitFunc->SetParNames("constant","MPV","sigma","2-Mip weight","3-Mip weight");
	    fitFunc->SetParameters(10,0.8,0.1,0.05,0.01);
	    fitFunc->SetParLimits(1,0.6,1.2);
	    fitFunc->SetParLimits(3,0,1);
	    fitFunc->SetParLimits(4,0,1);
	    
	  }
	  
	  
	  hEdist->Fit(fitFunc,"","",hEdist->GetBinCenter(hEdist->GetMaximumBin())-0.2,3);
	  fitFunc->Write(Form("FMD%d%c_etabin%d_fitfunc",det,ringChar,nEta),TObject::kWriteDelete);
	  
	}
	
	TH2F* hBg = pars->GetBackgroundCorrection(det,ringChar,0);
	EnergyDist->SetEnergyDistribution(det,ringChar,hBg->GetXaxis()->GetBinCenter(nEta),hEdist);
      }
      
    }
    
  }
  
  fin.Write();
  fin.Close();
  
  if(store) {
    TFile fcalib(pars->GetPath(pars->GetEdistID() ),"RECREATE");
    EnergyDist->Write(AliFMDAnaParameters::GetEdistID());
    fcalib.Close();
  }



}

//____________________________________________________________________
//
// EOF
//
