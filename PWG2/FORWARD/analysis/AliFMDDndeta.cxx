// Code to analyse dN/deta from the forward analysis
// This can plot the results 
// Also works for MC data 
//
// -- Author: Hans Hjersing Dalsgaard <canute@gmail.com>
#include "AliFMDDndeta.h"
#include "TFile.h"
#include "AliLog.h"
#include "TH1.h"
#include "AliFMDAnaParameters.h"
#include "AliFMDAnaCalibSharingEfficiency.h"
#include "TStyle.h"
//#include "TObjArray.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "iostream"
#include "TH3.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TF1.h"

#define SMALLNUMBER 0.0001

ClassImp(AliFMDDndeta)
//_____________________________________________________________________

AliFMDDndeta::AliFMDDndeta() 
: TObject(),
  fList(0),
  fMultList(),
  fNbinsToCut(2),
  fVtxCut1(-10),
  fVtxCut2(10),
  fIsInit(kFALSE),
  fIsGenerated(),
  fPrimEvents(),
  fEvents(),
  fPrimdNdeta(),
  fDrawAll(kFALSE)
{
  //AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  /* fDataObject = new TProfile3D("dataObject","dataObject",
			       pars->GetNetaBins(),-6,6,
			       20,0,2*TMath::Pi(),
			       pars->GetNvtxBins(),-0.5,pars->GetNvtxBins()-0.5);
  fDataObject->SetXTitle("#eta");
  fDataObject->SetYTitle("#varphi [radians]");
  fDataObject->SetZTitle("v_{z} [cm]");*/
  
  fAnalysisNames[0] = "Hits";
  fAnalysisNames[1] = "HitsTrVtx";
  fAnalysisNames[2] = "dNdeta";
  fAnalysisNames[3] = "dNdetaTrVtx";
  fAnalysisNames[4] = "dNdetaNSD";
  
  for(Int_t i=0; i<5;i++) 
    fMultList[i] = new TList();
}
//_____________________________________________________________________
void AliFMDDndeta::SetNames(Analysis what) {
  // Set names of histograms from analysis
  
  switch(what) {
  case kHits :
    fPrimEvents.Form("nMCEventsNoCuts"); //was nMCEvents
    fEvents.Form("nEvents");
    fPrimdNdeta.Form("hMultvsEtaNoCuts");
    break;
  case kHitsTrVtx :
    fPrimEvents.Form("nMCEvents"); 
    fEvents.Form("nEvents");
    fPrimdNdeta.Form("hMultvsEtaNoCuts");
    break;
  case kMult :
    fPrimEvents.Form("nMCEventsNoCuts"); 
    fEvents.Form("nEvents");
    fPrimdNdeta.Form("hMultvsEtaNoCuts");
    break;
  case kMultTrVtx :
    fPrimEvents.Form("nMCEvents"); 
    fEvents.Form("nEvents");
    fPrimdNdeta.Form("hMultvsEta");
    break;
  case kMultNSD :
    fPrimEvents.Form("nMCEventsNSDNoCuts"); 
    fEvents.Form("nNSDEvents");
    fPrimdNdeta.Form("hMultvsEtaNSDNoCuts");
    break;
    
  default:
    break;
  }
}
//_____________________________________________________________________
const char* AliFMDDndeta::GetAnalysisName(Analysis what, UShort_t det, Char_t ring, Int_t vtxbin) {
  //Get names of histograms
  const char* name = "";
  
  switch(what) {
  case kHits :
    name = Form("hits_NoCuts_FMD%d%c_vtxbin%d_proj",det,ring,vtxbin); //NoCuts added
    break;
  case kHitsTrVtx :
    name = Form("hits_FMD%d%c_vtxbin%d_proj",det,ring,vtxbin); 
    break;
  case kMult :
    name = Form("dNdeta_FMD%d%c_vtxbin%d_proj",det,ring,vtxbin);
    break;
  case kMultTrVtx :
    name = Form("dNdeta_FMD%d%c_TrVtx_vtxbin%d_proj",det,ring,vtxbin);
    break;
  case kMultNSD :
    name = Form("dNdetaNSD_FMD%d%c_vtxbin%d_proj",det,ring,vtxbin);
    break;
    
  default:
    break;
  }
  //std::cout<<name.Data()<<std::endl;
  return name;
}
//_____________________________________________________________________
const char* AliFMDDndeta::GetPrimName(Analysis what, UShort_t det, Char_t ring, Int_t vtxbin) {
  //Get names of primaries
  const char* name = "";
  
  switch(what) {
  case kHits :
    name = Form("hMCHits_nocuts_FMD%d%c_vtxbin%d",det,ring,vtxbin); //nocuts added
    break;
  case kHitsTrVtx :
    name = Form("hMCHits_FMD%d%c_vtxbin%d",det,ring,vtxbin);
    break;
  case kMult :
    name = Form("primmult_vtxbin%d",vtxbin);
    break;
  case kMultTrVtx :
    name = Form("primmult_NoCuts_vtxbin%d",vtxbin);
    break;
  case kMultNSD :
    name = Form("primmult_NSD_vtxbin%d",vtxbin);
    break;
  default:
    break;
  }
  return name;
}
//_____________________________________________________________________
void AliFMDDndeta::GenerateHits(Analysis what) {
  //Generate the hits distributions
  AliFMDAnaParameters* pars =  AliFMDAnaParameters::Instance();
  TH2F* hTmp         =  pars->GetBackgroundCorrection(1,'I',0);
  Int_t nVertexBins  =  pars->GetNvtxBins();
  Int_t nEtaBins     =  hTmp->GetNbinsX();
  Float_t etaMin     =  hTmp->GetXaxis()->GetXmin(); 
  Float_t etaMax     =  hTmp->GetXaxis()->GetXmax();

  for(Int_t i = 0; i<nVertexBins;i++) { // 
    TH1F* hHits = new TH1F(Form("hMCHits_vtxbin%d_%s",i,fAnalysisNames[what].Data()),Form("hHits_vtxbin%d_%s",i,fAnalysisNames[what].Data()),nEtaBins,etaMin,etaMax);
    hHits->Sumw2();
    fMultList[what]->Add(hHits);
  }
  
  for(Int_t det = 1; det<=3; det++) {
    Int_t maxRing = (det == 1 ? 0 : 1);
    for(Int_t ring = 0;ring<=maxRing;ring++) {
      
      Char_t ringChar = (ring == 0 ? 'I' : 'O');
      for(Int_t v=0; v< nVertexBins; v++) {
	TH1F* hits = (TH1F*)fList->FindObject(GetPrimName(what,det,ringChar,v));
	
	Int_t nNonZero = 0, nNonZeroInData = 0;
	
	//removing edges
	for(Int_t i =1;i<=hits->GetNbinsX();i++) {
	  if(hits->GetBinContent(i) !=0)
	    nNonZero++;
	}

	TH1F* sumMultHist = (TH1F*)fMultList[what]->FindObject(Form("hMCHits_vtxbin%d_%s",v,fAnalysisNames[what].Data()));
	for(Int_t i =1;i<=sumMultHist->GetNbinsX();i++) {
	  if(hits->GetBinContent(i) == 0 ) continue;
	  
	  nNonZeroInData++;
	  
	  if(nNonZeroInData <=fNbinsToCut || nNonZeroInData > (nNonZero - fNbinsToCut)) {
	    continue;
	  }
	  if(sumMultHist->GetBinContent(i) < SMALLNUMBER && TMath::Abs(hits->GetBinContent(i)) > 0){
	      sumMultHist->SetBinContent(i,hits->GetBinContent(i));
	      sumMultHist->SetBinError(i,hits->GetBinError(i));
	    }
	  else {
	      
	      
	    Float_t sumofw = (1/TMath::Power(hits->GetBinError(i),2)) + (1/TMath::Power(sumMultHist->GetBinError(i),2));
	    Float_t wav =  (hits->GetBinContent(i)*(1/TMath::Power(hits->GetBinError(i),2)) + sumMultHist->GetBinContent(i)*(1/TMath::Power(sumMultHist->GetBinError(i),2)))/sumofw;
	    Float_t error = 1/TMath::Sqrt(sumofw);
	    sumMultHist->SetBinContent(i, wav);
	    sumMultHist->SetBinError(i, error);
	  }
	}
      }
    }
  }
}

//_____________________________________________________________________
void AliFMDDndeta::Init(const Char_t* filename) {
  //Initialize everything from file
  TFile* fin = TFile::Open(filename);
  
  if(!fin) {
    AliWarning("No file - exiting !");
    return;
  }
  AliFMDAnaParameters* pars =  AliFMDAnaParameters::Instance();
  //pars->Init();
  
  TList* list = (TList*)fin->Get(Form("%s/BackgroundCorrected",pars->GetDndetaAnalysisName()));
  
  if(!list) //an old file ? Perhaps...
    list = (TList*)fin->Get("BackgroundCorrected");
  
  Init(list);
  
}
//_____________________________________________________________________
void AliFMDDndeta::Init(TList* list) {
  //Initialize everything from TList
    
  if(!list) {
    AliWarning("No list - exiting !");
    return;
  }
  
  fList = (TList*)list->Clone("inputList");
    
  fIsGenerated[kHits]      = kFALSE;
  fIsGenerated[kMult]      = kFALSE;  
  fIsGenerated[kMultTrVtx] = kFALSE;
  fIsGenerated[kHitsTrVtx] = kFALSE;
  fIsGenerated[kMultNSD]   = kFALSE;
  
  fIsInit = kTRUE;
}
//_____________________________________________________________________
void AliFMDDndeta::GenerateMult(Analysis what) {
  //Generate dNdeta
  
  if(!fIsInit) {
    AliWarning("Not initialised - call Init to remedy");
    return;
  }
  
  if(fIsGenerated[what]) {
    AliWarning("Already generated - have a look at the results!");
    return;
  }
  else fIsGenerated[what] = kTRUE;
  
  SetNames(what);
  
  if(what == kHits || what == kHitsTrVtx)
    GenerateHits(what);
  
  TH1I* hMCEvents = (TH1I*)fList->FindObject(fPrimEvents.Data());
  
  AliFMDAnaParameters* pars =  AliFMDAnaParameters::Instance();
  
  TH2F* hTmp         =  pars->GetBackgroundCorrection(1,'I',0);
  Int_t nVertexBins  =  pars->GetNvtxBins();
  Int_t nEtaBins     =  hTmp->GetNbinsX();
  Float_t etaMin     =  hTmp->GetXaxis()->GetXmin(); 
  Float_t etaMax     =  hTmp->GetXaxis()->GetXmax();
  
  for(Int_t i = 0; i<nVertexBins;i++) {
    TH1F* hMult = new TH1F(Form("hMult_vtxbin%d_%s",i,fAnalysisNames[what].Data()),Form("hMult_vtxbin%d_%s",i,fAnalysisNames[what].Data()),nEtaBins,etaMin,etaMax);
    hMult->Sumw2();
    fMultList[what]->Add(hMult);
  }
  
  for(Int_t det = 1; det<=3;det++) {
    Int_t maxRing = (det == 1 ? 0 : 1);
    for(Int_t iring = 0; iring<=maxRing; iring++) {
      Char_t ringChar = (iring == 0 ? 'I' : 'O');
      //Int_t nsec  = (iring == 0 ? 20 : 40);
      TH1F* hRingMult= new TH1F(Form("hRingMult_FMD%d%c_%s",det,ringChar,fAnalysisNames[what].Data()),Form("hRingMult_FMD%d%c_%s",det,ringChar,fAnalysisNames[what].Data()),nEtaBins,etaMin,etaMax);
      fMultList[what]->Add(hRingMult);
      //      TProfile* phiprofile     = new TProfile(Form("dNdphiFMD%d%c",det,ringChar), Form("dNdphiFMD%d%c;#Phi",det,ringChar), nsec , 0, 2*TMath::Pi());
      //fMultList[what]->Add(phiprofile);
    }
  }
  TH1I*  hEvents  = (TH1I*)fList->FindObject(fEvents.Data());
  TH1F*  hPrimVtx = 0;
  
  	
    
  for(Int_t det = 1; det<=3; det++) {
    Int_t maxRing = (det == 1 ? 0 : 1);
    for(Int_t ring = 0;ring<=maxRing;ring++) {
      
      Char_t ringChar = (ring == 0 ? 'I' : 'O');
      for(Int_t v=0; v< nVertexBins; v++) {
	if(det == 1) {
	  if(what == kHits || what == kHitsTrVtx)
	    hPrimVtx = (TH1F*)fMultList[what]->FindObject(Form("hMCHits_vtxbin%d_%s",v,fAnalysisNames[what].Data()));
	  else
	    hPrimVtx = (TH1F*)fList->FindObject(GetPrimName(what,det,ringChar,v));
	  

	  Float_t nPrimevents = hMCEvents->GetBinContent(v+1);
	  Float_t xb = hPrimVtx->GetNbinsX();
	  Float_t xr = hPrimVtx->GetXaxis()->GetXmax() - hPrimVtx->GetXaxis()->GetXmin(); 
	  hPrimVtx->Scale(xb / xr );
	  if(nPrimevents > 0)
	    hPrimVtx->Scale(1/nPrimevents);
	  
	}
	
	Double_t delta = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
	Float_t vtxZ1   = (delta*v) - pars->GetVtxCutZ();
	Float_t vtxZ2   = (delta*(v+1)) - pars->GetVtxCutZ();
	
	if(vtxZ1< fVtxCut1 || vtxZ2 >fVtxCut2)
	  continue;
	Float_t nEvents = (Float_t)hEvents->GetBinContent(v+1);
	
	//TH1F* multhistproj = (TH1F*)fList->FindObject(Form("dNdeta_FMD%d%c_vtxbin%d_proj",det,ringChar,v));
       	

	
	
	
	
	TH1F* multhistproj = (TH1F*)fList->FindObject(GetAnalysisName(what,det,ringChar,v));
	if(nEvents)
	  multhistproj->Scale(1/nEvents);
	Float_t xnBins = multhistproj->GetNbinsX();
	Float_t xrange = multhistproj->GetXaxis()->GetXmax() - multhistproj->GetXaxis()->GetXmin(); 
	multhistproj->Scale(xnBins / xrange);
	
	//TH2F* multhist = (TH2F*)fList->FindObject(Form("dNdeta_FMD%d%c_vtxbin%d",det,ringChar,v));

	
	//if(nEvents)
	//  multhist->Scale(1/nEvents);
	/*	for(Int_t i=1;i<=multhist->GetNbinsX();i++) {
	  for(Int_t j=1;j<=multhist->GetNbinsY();j++) {
	    
	    if (multhist->GetBinContent(i,j) <= 0.0001) continue;
	    //std::cout<<multhist->GetBinContent(i,j)<<std::endl;
	    fDataObject->Fill(multhist->GetXaxis()->GetBinCenter(i),
			      multhist->GetYaxis()->GetBinCenter(j),
			      v,
			      multhist->GetBinContent(i,j), multhist->GetBinError(i,j));
	    //fDataObject->SetBinContent(i,j,v,multhist->GetBinContent(i,j));
	    //fDataObject->SetBinError(i,j,v,multhist->GetBinError(i,j));
	  } 
	}
	*/	
	//if(nEvents)
	//  multhist->Scale(1/nEvents);
	//if(ringChar == 'O')
	//  multhist->RebinY(2);
	//removing edges
	
	Int_t nNonZero = 0, nNonZeroInData = 0;
	
	for(Int_t i =1;i<=multhistproj->GetNbinsX();i++) {
	  if(multhistproj->GetBinContent(i) !=0)
	    nNonZero++;
	}
	Int_t nBinsOld = fNbinsToCut;
	//if(det == 1 && ringChar =='I') {
	//  fNbinsToCut = 0;
	//	}
	TH1F* hRingMult = (TH1F*)fMultList[what]->FindObject(Form("hRingMult_FMD%d%c_%s",det,ringChar,fAnalysisNames[what].Data()));
	
	for(Int_t i=1; i<=hRingMult->GetNbinsX(); i++) {
	  if(multhistproj->GetBinContent(i)!=0) {
	    nNonZeroInData++;
	    if(nNonZeroInData <=fNbinsToCut || nNonZeroInData > (nNonZero - fNbinsToCut)) {
	      continue;
	    }
	    Float_t oldweight = 0;
	    Float_t oldwav    = 0;
	    if(hRingMult->GetBinError(i)>0) {
	      oldweight = 1/TMath::Power(hRingMult->GetBinError(i),2);
	      oldwav    = oldweight*hRingMult->GetBinContent(i);
	    }
	    Float_t  weight = 1/TMath::Power(multhistproj->GetBinError(i),2);
	    Float_t  wav    = oldwav + weight*multhistproj->GetBinContent(i);
	    Float_t  sumofw = oldweight + weight;
	    if(sumofw) {
	      Float_t error = 1/TMath::Sqrt(sumofw);
	      
	      hRingMult->SetBinContent(i,wav/sumofw);
	      hRingMult->SetBinError(i,error);
	      
	    }
	  }
	}
	nNonZeroInData = 0;
	
	TH1F* sumMultHist = (TH1F*)fMultList[what]->FindObject(Form("hMult_vtxbin%d_%s",v,fAnalysisNames[what].Data()));
	
	for(Int_t i =1;i<=sumMultHist->GetNbinsX();i++) {
	  if(multhistproj->GetBinContent(i) != 0 ) {	  
	    nNonZeroInData++;
	  
	    if(nNonZeroInData <=fNbinsToCut || nNonZeroInData > (nNonZero - fNbinsToCut)) {
	      continue;
	    }
	    if(sumMultHist->GetBinContent(i) < SMALLNUMBER && TMath::Abs(multhistproj->GetBinContent(i)) > 0){
	      sumMultHist->SetBinContent(i,multhistproj->GetBinContent(i));
	      sumMultHist->SetBinError(i,multhistproj->GetBinError(i));
	    }
	    else {
	      
	      
	      Float_t sumofw = (1/TMath::Power(multhistproj->GetBinError(i),2)) + (1/TMath::Power(sumMultHist->GetBinError(i),2));
	      Float_t wav =  (multhistproj->GetBinContent(i)*(1/TMath::Power(multhistproj->GetBinError(i),2)) + sumMultHist->GetBinContent(i)*(1/TMath::Power(sumMultHist->GetBinError(i),2)))/sumofw;
	      Float_t error = 1/TMath::Sqrt(sumofw);
	      sumMultHist->SetBinContent(i, wav);
	      sumMultHist->SetBinError(i, error);
	    }
	    
	
	  }
	}
	fNbinsToCut = nBinsOld;
      }
      
    }
  }
  
  TH1F* sumMult  = new TH1F(Form("dNdeta_%s",fAnalysisNames[what].Data()),Form("dNdeta_%s",fAnalysisNames[what].Data()),nEtaBins,etaMin,etaMax);
  sumMult->Sumw2();
  fMultList[what]->Add(sumMult);
  TH1F* primMult  = new TH1F(Form("primary_%s",fAnalysisNames[what].Data()),Form("primary_%s",fAnalysisNames[what].Data()),nEtaBins,etaMin,etaMax);
  primMult->Sumw2();
  fMultList[what]->Add(primMult);

  Float_t wav  = 0, sumofw=0, weight = 0, primwav  = 0, primsumofw=0, primweight = 0;//, etaofbin = 0;
  for(Int_t i =1; i<=sumMult->GetNbinsX();i++) {
    
    for(Int_t v = 0; v<nVertexBins;v++) {
      if(what == kHits || what == kHitsTrVtx)
	hPrimVtx = (TH1F*)fMultList[what]->FindObject(Form("hMCHits_vtxbin%d_%s",v,fAnalysisNames[what].Data()));
      else
	hPrimVtx = (TH1F*)fList->FindObject(GetPrimName(what,0,0,v));
      TH1F* sumMultHist = (TH1F*)fMultList[what]->FindObject(Form("hMult_vtxbin%d_%s",v,fAnalysisNames[what].Data()));
      //etaofbin += sumMultHist->GetBinCenter(i);
      // if( hPrimVtx->GetBinContent(i)!=0 && hPrimVtx->GetBinError(i)>0.001 && sumMultHist->GetBinContent(i)!=0) {
      if( TMath::Abs(hPrimVtx->GetBinContent(i)) > 0) {
	primweight = 1/TMath::Power(hPrimVtx->GetBinError(i),2);
	primwav    = primwav + primweight*hPrimVtx->GetBinContent(i);
	primsumofw = primsumofw + primweight;
      }
      if(TMath::Abs(sumMultHist->GetBinContent(i)) > 0) {
	
	weight = 1/TMath::Power(sumMultHist->GetBinError(i),2);
	wav    = wav + weight*sumMultHist->GetBinContent(i);
	sumofw = sumofw + weight;
      }
      
    }
    if( primsumofw !=0 ) {// && sumofw !=0) {
      Float_t primerror = 1/TMath::Sqrt(primsumofw);
      
      primMult->SetBinContent(i,primwav/primsumofw);
      primMult->SetBinError(i,primerror);
    }
    else {
      primMult->SetBinContent(i,0);
      primMult->SetBinError(i,0);
    }
    
    if(sumofw) {
      
      Float_t error = 1/TMath::Sqrt(sumofw);
      
      sumMult->SetBinContent(i,wav/sumofw);
      sumMult->SetBinError(i,error);
    }
    
    
    wav = 0;
    sumofw = 0;
    weight = 0;
    primwav = 0;
    primsumofw = 0;
    primweight = 0;
    
  }
  
}
//_____________________________________________________________________
void AliFMDDndeta::DrawDndeta(Analysis what, Int_t rebin, Bool_t realdata, TString pythiafile) {
  //Draw everything
  if(!fIsInit) {
    AliWarning("Not initialised - call Init to remedy");
    return;
  }

  if(!fIsGenerated[what]) {
    AliWarning("Not generated - generate first then draw!");
    return;
  }
  AliFMDAnaParameters* pars =  AliFMDAnaParameters::Instance();
  SetNames(what);
  //UA5 data NSD
  Float_t x[19] = {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,
		   2.625,2.875,3.125,3.375,3.625,3.875,4.125,4.375,4.625};
  Float_t x2[19] = {-0.125,-0.375,-0.625,-0.875,-1.125,-1.375,-1.625,-1.875,
		    -2.125,-2.375,-2.625,-2.875,-3.125,-3.375,-3.625,-3.875,
		    -4.125,-4.375,-4.625};
  
  
  Float_t y[19] = {3.48, 3.38, 3.52, 3.68, 3.71, 3.86, 3.76, 3.66, 3.72 ,3.69,
		   3.56, 3.41, 3.15, 3.09, 2.74, 2.73, 2.32, 1.99, 1.69};
  
  Float_t ey[19] = {0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,
		    0.07,0.07,0.08,0.08,0.09,0.09,0.1,0.13};
    
  TGraphErrors* graph = new TGraphErrors(19,x,y,NULL,ey);
  TGraphErrors* graph2 = new TGraphErrors(19,x2,y,NULL,ey);
  //UA5 data INEL
  Float_t xinel[19] = {0.125, 0.375, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.125, 4.375, 4.625};
  Float_t xinel2[19] = {-0.125, -0.375, -0.625, -0.875, -1.125, -1.375, -1.625, -1.875, -2.125, -2.375, -2.625, -2.875, -3.125, -3.375, -3.625, -3.875, -4.125, -4.375, -4.625};
  Float_t yinel[19] = {3.14, 3.04, 3.17, 3.33, 3.33, 3.53, 3.46, 3.41, 3.45, 3.39, 3.07, 3.07, 2.93, 2.93, 2.55, 2.48, 2.18, 1.91, 1.52};
  Float_t eyinel[19] = {0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.08, 0.08, 0.09, 0.09, 0.1, 0.13};
 TGraphErrors* graphinel = new TGraphErrors(19,xinel,yinel,NULL,eyinel);
 TGraphErrors* graphinel2 = new TGraphErrors(19,xinel2,yinel,NULL,eyinel);
  
  graph->SetMarkerStyle(22);
  graph->SetMarkerColor(kBlue);
  graphinel->SetMarkerStyle(20);
  graphinel->SetMarkerColor(kGreen);
  graphinel2->SetMarkerStyle(24);
  graphinel2->SetMarkerColor(kGreen);
  //graph->Draw("P");
  graph2->SetMarkerStyle(26);
  graph2->SetMarkerColor(kBlue);
  graphinel->GetHistogram()->SetStats(kFALSE);
  graphinel2->GetHistogram()->SetStats(kFALSE);
  
  
  //End UA5 data
  
  //Published ALICE data
  // Plot: p7742_d1x1y1
  
  // INEL
  
  TGraphAsymmErrors* p7742D1x1y1 = 0;
  double p7742D1x1y1xval[] = { -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 
    0.5, 0.7, 0.9, 1.1, 1.3 };
  double p7742D1x1y1xerrminus[] = { 0.09999999999999987, 0.09999999999999987, 0.09999999999999998, 0.10000000000000009, 0.09999999999999998, 0.10000000000000003, 0.1, 0.1, 0.09999999999999998, 
    0.09999999999999998, 0.09999999999999998, 0.09999999999999998, 0.10000000000000009, 0.10000000000000009 };
  double p7742D1x1y1xerrplus[] = { 0.10000000000000009, 0.10000000000000009, 0.09999999999999998, 0.09999999999999998, 0.09999999999999998, 0.09999999999999998, 0.1, 0.1, 0.10000000000000003, 
    0.09999999999999998, 0.10000000000000009, 0.09999999999999998, 0.09999999999999987, 0.09999999999999987 };
  double p7742D1x1y1yval[] = { 3.28, 3.28, 3.22, 3.12, 3.06, 3.02, 2.98, 3.02, 3.02, 
    3.05, 3.15, 3.21, 3.26, 3.33 };
  double p7742D1x1y1yerrminus[] = { 0.06324555320336758, 0.06324555320336758, 0.06324555320336758, 0.06324555320336758, 0.06324555320336758, 0.05385164807134505, 0.05385164807134505, 0.05385164807134505, 0.05385164807134505, 
    0.06324555320336758, 0.06324555320336758, 0.06324555320336758, 0.06324555320336758, 0.06324555320336758 };
  double p7742D1x1y1yerrplus[] = { 0.08246211251235322, 0.08246211251235322, 0.08246211251235322, 0.08246211251235322, 0.08246211251235322, 0.08246211251235322, 0.07280109889280519, 0.08246211251235322, 0.08246211251235322, 
    0.08246211251235322, 0.08246211251235322, 0.08246211251235322, 0.08246211251235322, 0.08246211251235322 };
  int p7742D1x1y1numpoints = 14;
  p7742D1x1y1 = new TGraphAsymmErrors(p7742D1x1y1numpoints, p7742D1x1y1xval, p7742D1x1y1yval, p7742D1x1y1xerrminus, p7742D1x1y1xerrplus, p7742D1x1y1yerrminus, p7742D1x1y1yerrplus);
  p7742D1x1y1->SetName("/HepData/7742/d1x1y1");
  p7742D1x1y1->SetTitle("/HepData/7742/d1x1y1");
  p7742D1x1y1->SetMarkerStyle(21);
  p7742D1x1y1->SetMarkerColor(kRed);
  // p7742D1x1y1->Draw("Psame");
  
  
  //NSD
  TGraphAsymmErrors*  p7742D2x1y1 = 0;
 double p7742D2x1y1xval[] = { -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 
    0.5, 0.7, 0.9, 1.1, 1.3 };
  double p7742D2x1y1xerrminus[] = { 0.09999999999999987, 0.09999999999999987, 0.09999999999999998, 0.10000000000000009, 0.09999999999999998, 0.10000000000000003, 0.1, 0.1, 0.09999999999999998, 
    0.09999999999999998, 0.09999999999999998, 0.09999999999999998, 0.10000000000000009, 0.10000000000000009 };
  double p7742D2x1y1xerrplus[] = { 0.10000000000000009, 0.10000000000000009, 0.09999999999999998, 0.09999999999999998, 0.09999999999999998, 0.09999999999999998, 0.1, 0.1, 0.10000000000000003, 
    0.09999999999999998, 0.10000000000000009, 0.09999999999999998, 0.09999999999999987, 0.09999999999999987 };
  double p7742D2x1y1yval[] = { 3.9, 3.89, 3.81, 3.7, 3.64, 3.59, 3.53, 3.58, 3.59, 
    3.61, 3.74, 3.8, 3.87, 3.95 };
  double p7742D2x1y1yerrminus[] = { 0.13341664064126335, 0.13152946437965907, 0.13152946437965907, 0.1216552506059644, 0.1216552506059644, 0.1216552506059644, 0.1216552506059644, 0.1216552506059644, 0.1216552506059644, 
    0.1216552506059644, 0.1216552506059644, 0.13152946437965907, 0.13152946437965907, 0.13341664064126335 };
  double p7742D2x1y1yerrplus[] = { 0.13341664064126335, 0.13152946437965907, 0.13152946437965907, 0.1216552506059644, 0.1216552506059644, 0.1216552506059644, 0.1216552506059644, 0.1216552506059644, 0.1216552506059644, 
    0.1216552506059644, 0.1216552506059644, 0.13152946437965907, 0.13152946437965907, 0.13341664064126335 };
  int p7742D2x1y1numpoints = 14;

  p7742D2x1y1 = new TGraphAsymmErrors(p7742D2x1y1numpoints, p7742D2x1y1xval, p7742D2x1y1yval, p7742D2x1y1xerrminus, p7742D2x1y1xerrplus, p7742D2x1y1yerrminus, p7742D2x1y1yerrplus);
  p7742D2x1y1->SetName("/HepData/7742/d2x1y1");
  p7742D2x1y1->SetTitle("/HepData/7742/d2x1y1");
  p7742D2x1y1->SetMarkerStyle(21);
  p7742D2x1y1->SetMarkerColor(kRed);
  //p7742D2x1y1.Draw("AP");






  //End official data
  
  // CMS published NSD data
  
  TGraphAsymmErrors* p7743D8x1y1 = 0;
  double p7743D8x1y1xval[] = { -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 
    2.25 };
  double p7743D8x1y1xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25 };
  double p7743D8x1y1xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25 };
  double p7743D8x1y1yval[] = { 3.6, 3.73, 3.62, 3.54, 3.48, 3.48, 3.54, 3.62, 3.73,  3.6 };
  double p7743D8x1y1yerrminus[] = { 0.13, 0.14, 0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 0.14,0.13 };
  double p7743D8x1y1yerrplus[] = { 0.13, 0.14, 0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 0.14, 0.13 };
  int p7743D8x1y1numpoints = 10;
  p7743D8x1y1 = new TGraphAsymmErrors(p7743D8x1y1numpoints, p7743D8x1y1xval, p7743D8x1y1yval, p7743D8x1y1xerrminus, p7743D8x1y1xerrplus, p7743D8x1y1yerrminus, p7743D8x1y1yerrplus);

  p7743D8x1y1->SetMarkerStyle(20);
  p7743D8x1y1->SetMarkerColor(kBlack);
  
  //End CMS
  
  
  TH1I* hEvents   = (TH1I*)fList->FindObject(fEvents.Data());
  
  TH1I* hMCEvents = (TH1I*)fList->FindObject(fPrimEvents.Data());
  
   //SPD part
  TH1D* hSPDana = (TH1D*)fList->FindObject(Form("dNdeta_SPD_vtxbin%d_proj",5));
  
  for(Int_t nn = 0; nn < pars->GetNvtxBins() ; nn++) {
    TH1D* hSPDanalysis = (TH1D*)fList->FindObject(Form("dNdeta_SPD_vtxbin%d_proj",nn));
    TH1D* hSPDanalysisTrVtx = (TH1D*)fList->FindObject(Form("dNdetaTrVtx_SPD_vtxbin%d_proj",nn));
    TH1D* hSPDanalysisNSD = (TH1D*)fList->FindObject(Form("dNdetaNSD_SPD_vtxbin%d_proj",nn));
    
    Float_t nEventsSPD = (Float_t)hEvents->GetBinContent(nn+1);
    if(!nEventsSPD) continue; 
    hSPDanalysis->Scale(1/nEventsSPD);
    hSPDanalysisTrVtx->Scale(1/nEventsSPD);
    hSPDanalysisNSD->Scale(1/nEventsSPD);
  }
  Float_t lowlimits[10] = {-1,-1.25,-1.5,-1.75,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9};
  Float_t highlimits[10] = {1.9,1.9,1.9,1.9,1.9,1.9,1.75,1.5,1.25,1};

  TH1F* hSPDcombi = new TH1F("SPDcombi","SPDcombi",hSPDana->GetNbinsX(),hSPDana->GetXaxis()->GetXmin(),hSPDana->GetXaxis()->GetXmax());
  TH1D* hSPDanalysis = 0;
  for(Int_t kk = 1; kk <=hSPDana->GetNbinsX(); kk++) {
    
    Float_t weight = 0, wav=0,sumofw = 0;
    
    // if(hSPDana->GetBinCenter(kk) < lowlimits[nn])
    //	continue;
    //if(hSPDana->GetBinCenter(kk) > highlimits[nn])
    //  continue;
    
    for(Int_t nn =0; nn < pars->GetNvtxBins() ; nn++) {   
      Double_t delta = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
      Float_t vtxZ1   = (delta*nn) - pars->GetVtxCutZ();
      Float_t vtxZ2   = (delta*(nn+1)) - pars->GetVtxCutZ();
      
      if(vtxZ1<fVtxCut1 || vtxZ2 >fVtxCut2)
	continue;
      
      //std::cout<<vtxZ1<<"   "<<vtxZ2<<std::endl;
      
      if(what == kMult)
	hSPDanalysis = (TH1D*)fList->FindObject(Form("dNdeta_SPD_vtxbin%d_proj",nn));
      else if(what == kMultNSD)
	hSPDanalysis = (TH1D*)fList->FindObject(Form("dNdetaNSD_SPD_vtxbin%d_proj",nn));
      else if(what == kMultTrVtx)
	hSPDanalysis = (TH1D*)fList->FindObject(Form("dNdetaTrVtx_SPD_vtxbin%d_proj",nn));
      else continue;
      if(hSPDanalysis->GetBinCenter(kk) < lowlimits[nn])
    	continue;
      if(hSPDanalysis->GetBinCenter(kk) > highlimits[nn])
	continue;
      //std::cout<<hSPDanalysis->GetBinCenter(kk)<<"  "<<lowlimits[nn]<<std::endl;
     
      Float_t mult = hSPDanalysis->GetBinContent(kk);
      Float_t error = hSPDanalysis->GetBinError(kk);
      /*
      if(mult > 0 && hSPDanalysis->GetBinContent(kk-1) < SMALLNUMBER)
	continue;
      if(mult > 0 && hSPDanalysis->GetBinContent(kk-2) < SMALLNUMBER)
	continue;
      
      if(mult > 0 && hSPDanalysis->GetBinContent(kk+1) < SMALLNUMBER)
	continue;
      if(mult > 0 && hSPDanalysis->GetBinContent(kk+2) < SMALLNUMBER)
	continue;
      */
      if(mult > 0) {
	weight = 1/TMath::Power(error,2);
	wav    = wav + weight*mult;
	sumofw = sumofw + weight;
      }
      
    }
    
    if(sumofw && wav) {
      Float_t errorTotal = 1/TMath::Sqrt(sumofw);
      hSPDcombi->SetBinContent(kk,wav/sumofw);
      hSPDcombi->SetBinError(kk,errorTotal);
    }
  }
    
  
  Float_t xb1 = hSPDcombi->GetNbinsX();
  Float_t xr1 = hSPDcombi->GetXaxis()->GetXmax() - hSPDcombi->GetXaxis()->GetXmin(); 
  hSPDcombi->Scale(xb1 / xr1);
  //hSPDcombi->Rebin(rebin);
  //hSPDcombi->Scale(1/(Float_t)rebin);
  //RebinHistogram(hSPDcombi,rebin);
  hSPDcombi->SetMarkerStyle(29);
  hSPDcombi->SetMarkerColor(kBlue);
  //hSPDcombi->DrawCopy("same");
  
  
  TCanvas* c1 = new TCanvas("cvertexbins","Overlaps",1400,800);
  c1->Divide(5,2);
  TCanvas* cCorrections = 0;
  TCanvas* cCorrectionsPhi = 0;
  if(fDrawAll) {
    cCorrections = new TCanvas("corrections","Correction vs Eta",1400,800);
    cCorrections->Divide(5,2);
    
    cCorrectionsPhi = new TCanvas("correctionsPhi","Correction vs Phi",1400,800);
    cCorrectionsPhi->Divide(5,2);
  }
  //TCanvas* cphi = new TCanvas("dNdphi","dNdphi",1280,1024);
  // cphi->Divide(3,2);
  
  Int_t nVertexBins  =  pars->GetNvtxBins();
  
  
  //Int_t npadphi = 0;
  TH1F*  hPrimVtx = 0;
  for(Int_t det = 1; det<=3; det++) {
    Int_t maxRing = (det == 1 ? 0 : 1);
    for(Int_t ring = 0;ring<=maxRing;ring++) {
      //  npadphi++;
      
      for(Int_t v=0; v< nVertexBins; v++) {
	Char_t ringChar = (ring == 0 ? 'I' : 'O');
	Double_t delta = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
	Float_t vtxZ1   = (delta*v) - pars->GetVtxCutZ();
	Float_t vtxZ2   = (delta*(v+1)) - pars->GetVtxCutZ();
	  
	if(vtxZ1<fVtxCut1 || vtxZ2 >fVtxCut2)
	  continue;
	Float_t nEvents = (Float_t)hEvents->GetBinContent(v+1);
	if(fDrawAll) {
	TH2F* hBgCor  = pars->GetBackgroundCorrection(det,ringChar,v);
	if(hBgCor)  {
	  TH1D* hBgProj = hBgCor->ProjectionX(Form("hBgProj_FMD%d%c_vtxbin%d",det,ringChar,v));
	  TH1D* hBgProjPhi = hBgCor->ProjectionY(Form("hBgProjPhi_FMD%d%c_vtxbin%d",det,ringChar,v));
	  
	  hBgProjPhi->SetName(Form("FMD%d%c",det,ringChar));
	  hBgProjPhi->SetTitle("");//Form("FMD%d",det));
	  hBgProj->SetTitle("");
	  //Float_t scalefactor = (hBgProj->GetXaxis()->GetXmax() - hBgProj->GetXaxis()->GetXmin()) / (Float_t)hBgProj->GetNbinsX();
	  Float_t scalefactor = 1/(Float_t)hBgCor->GetNbinsY();
	  
	  Float_t nFilledEtaBins = 0;
	  
	  for(Int_t jj = 1; jj<=hBgProj->GetNbinsX(); jj++) {
	    if(hBgProj->GetBinContent(jj) > 0.01)
	      nFilledEtaBins++;
	  }
	  
	  Float_t scalefactorPhi = 1/nFilledEtaBins;
	  
	  
	  hBgProj->Scale(scalefactor);
	  cCorrections->cd(v+1);
	  
	  hBgProj->SetMarkerColor(det + 2*ring);
	  hBgProj->SetMarkerStyle(19 + det + 2*ring );
	  
	  if((det + 2*ring) == 5)  hBgProj->SetMarkerColor(det + 2*ring+1);
	  if((19 + det + 2*ring ) == 24)  hBgProj->SetMarkerStyle(29 );
	  
	  if(det == 1) {
	    hBgProj->GetYaxis()->SetRangeUser(0,4);
	    hBgProj->SetStats(kFALSE);
	    hBgProj->SetXTitle("#eta");
	    hBgProj->DrawCopy("PE");
	    TLatex* l = new TLatex(0.14,0.92,Form("Vtx range [%.1f, %.1f]",vtxZ1,vtxZ2));
	    l->SetNDC(kTRUE);
	    l->Draw();
	  }
	  else
	    hBgProj->DrawCopy("PEsame");
	
	  hBgProjPhi->Scale(scalefactorPhi);
	  cCorrectionsPhi->cd(v+1);
	  hBgProjPhi->SetMarkerColor(det + 2*ring);
	  hBgProjPhi->SetMarkerStyle(19 + det + 2*ring );
	  if((det + 2*ring) == 5)  hBgProjPhi->SetMarkerColor(det + 2*ring+1);
	  if((19 + det + 2*ring ) == 24)  hBgProjPhi->SetMarkerStyle(29 );
	  if(det == 1) {
	    hBgProjPhi->GetYaxis()->SetRangeUser(1,5);
	    hBgProjPhi->SetStats(kFALSE);
	    
	    hBgProjPhi->SetXTitle("#Phi");
	    hBgProjPhi->DrawCopy("PE");
	   
	    
	  }
	  else
	    hBgProjPhi->DrawCopy("PEsame");
	  
	}
	
	}
	
	TH1F* multhistproj = (TH1F*)fList->FindObject(GetAnalysisName(what,det,ringChar,v));
		
	c1->cd(v+1);
	
	if(det==1)  {
	  multhistproj->SetMarkerStyle(20); 
	  multhistproj->SetMarkerColor(1); 
	}
	if(det==2 && ringChar=='I') {
	  multhistproj->SetMarkerStyle(21);
	  multhistproj->SetMarkerColor(2); 
	}
	if(det==2 && ringChar=='O') {
	  multhistproj->SetMarkerStyle(3);
	  multhistproj->SetMarkerColor(3); 
	}
	if(det==3 && ringChar=='I') {
	  multhistproj->SetMarkerStyle(22);
	  multhistproj->SetMarkerColor(4); 
	}
	if(det==3 && ringChar=='O') {
	  multhistproj->SetMarkerStyle(23);
	  multhistproj->SetMarkerColor(6); 
	}
	
	
	
	if(det == 1) {
	  if(what == kHits || what == kHitsTrVtx)
	    hPrimVtx = (TH1F*)fMultList[what]->FindObject(Form("hMCHits_vtxbin%d_%s",v,fAnalysisNames[what].Data()));
	  else
	    hPrimVtx = (TH1F*)fList->FindObject(GetPrimName(what,det,ringChar,v));
	  // std::cout<<hPrimVtx<<"   "<<kHits<<"   "<<kHitsTrVtx<<"   "<<what<<std::endl;
	  //std::cout<<Form("hMCHits_vtxbin%d_%s",v,fAnalysisNames[what].Data())<<std::endl;
	  hPrimVtx->SetTitle("");
	  TLatex* l = new TLatex(0.14,0.92,Form("Vtx range [%.1f, %.1f], %d events",vtxZ1,vtxZ2,(Int_t)nEvents));
	  l->SetNDC(kTRUE);
	  hPrimVtx->SetFillColor(kGray);
	  hPrimVtx->SetLabelFont(132,"xyz");
	  hPrimVtx->SetStats(0000);
	  hPrimVtx->GetXaxis()->SetTitle("#eta");
	  hPrimVtx->GetYaxis()->SetRangeUser(0,6);
	  hPrimVtx->DrawCopy("E3");
	  l->Draw();
	  multhistproj->DrawCopy("same");
	  TH1D* hSPDanavtxbin = 0;
	  if(what == kMult || what == kMultTrVtx || what == kMultNSD)  {
	  
	  if(what == kMult)
	    hSPDanavtxbin = (TH1D*)fList->FindObject(Form("dNdeta_SPD_vtxbin%d_proj",v));
	  if(what == kMultTrVtx)
	    hSPDanavtxbin = (TH1D*)fList->FindObject(Form("dNdetaTrVtx_SPD_vtxbin%d_proj",v));
	  
	  if(what == kMultNSD)
	    hSPDanavtxbin = (TH1D*)fList->FindObject(Form("dNdetaNSD_SPD_vtxbin%d_proj",v));
	  
	  hSPDanavtxbin->SetMarkerColor(kBlue);
	  hSPDanavtxbin->SetMarkerStyle(kStar);
	  hSPDanavtxbin->Scale(xb1 / xr1);
	  hSPDanavtxbin->DrawCopy("same");
	  
	  if(what != kMultNSD) {
	    graphinel->Draw("sameP");
	    graphinel2->Draw("sameP");
	  }
	  else
	   {
	     graph->Draw("sameP");
	     graph2->Draw("sameP");
	   }
	  }
	}
	else
	  multhistproj->DrawCopy("same");
	
      
      }
    }
  }
  
  //Legends for corrections
  if(fDrawAll) {
    for(Int_t v=0; v< nVertexBins; v++) {
      TPad* pad= (TPad*)cCorrectionsPhi->cd(v+1);
      pad->BuildLegend(0.15,0.45,0.45,0.9);
      Double_t delta = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
      Float_t vtxZ1   = (delta*v) - pars->GetVtxCutZ();
      Float_t vtxZ2   = (delta*(v+1)) - pars->GetVtxCutZ();
      TLatex* l = new TLatex(0.14,0.92,Form("Vtx range [%.1f, %.1f]",vtxZ1,vtxZ2));
      l->SetNDC(kTRUE);
      l->Draw();
    }
  }
  for(Int_t v=0; v< nVertexBins; v++) {
    TH1F* sumMultHist = (TH1F*)fMultList[what]->FindObject(Form("hMult_vtxbin%d_%s",v,fAnalysisNames[what].Data()));
    c1->cd(v+1);
    sumMultHist->SetMarkerStyle(25);
    sumMultHist->DrawCopy("same");
    
  }
  TH1F* primMult = (TH1F*)fMultList[what]->FindObject(Form("primary_%s",fAnalysisNames[what].Data()));
  TH1F* sumMult  = (TH1F*)fMultList[what]->FindObject(Form("dNdeta_%s",fAnalysisNames[what].Data()));
  sumMult->SetMarkerStyle(23);
  sumMult->SetMarkerColor(kRed);
  for(Int_t det = 1; det<=3;det++) {
    Int_t maxRing = (det == 1 ? 0 : 1);
    for(Int_t iring = 0; iring<=maxRing; iring++) {
      Char_t ringChar = (iring == 0 ? 'I' : 'O');
      TH1F* hRingMult= (TH1F*)fMultList[what]->FindObject(Form("hRingMult_FMD%d%c_%s",det,ringChar,fAnalysisNames[what].Data()));
      hRingMult->Add(primMult,-1);
      hRingMult->Divide(primMult);
    }
  }
  
  
  // End of SPD
  
  
  //TH1F*  hPrim = (TH1F*)fList->FindObject(fPrimEvents.Data());
  //  TFile testhit("/home/canute/ALICE/Simulations/TestOfAnalysis2/hitsdist.root");
  // TH1F* hHitCor = (TH1F*)testhit.Get("selectedHits");
  // hHitCor->Sumw2();
  // sumMult->Divide(hHitCor);
  TH1F*  hPrim = (TH1F*)fList->FindObject(fPrimdNdeta.Data());
  Float_t nPrimevents = hMCEvents->GetEntries();
  //std::cout<<hMCEvents->GetEntries()<<std::endl;
  Float_t xb = hPrim->GetNbinsX();
  Float_t xr = hPrim->GetXaxis()->GetXmax() - hPrim->GetXaxis()->GetXmin(); 
  //std::cout<<xb/xr<<std::endl;
  hPrim->Scale(xb / xr );
  if(nPrimevents > 0)
    hPrim->Scale(1/nPrimevents);
  
  
  //  primMult->Rebin(rebin);
  // primMult->Scale(1/rebin);
  // hPrim->Rebin(rebin);
  // hPrim->Scale(1/rebin);
  if(what != kHits && what != kHitsTrVtx)
    primMult = hPrim;
  
  //  hReldif->Add(hPrim,-1);
  // hReldif->Divide(hPrim);
  


  /////////////*******************New thing!!!
  
  //TH3D* hist3d = (TH3D*)fDataObject->ProjectionXYZ("hist3d");
  // fDataObject->Sumw2();
  //TH2F* projeta = (TH2F*)fDataObject->Project3D("yx");
  
  //  TProfile2D* projeta = fDataObject->Project3DProfile("yx");
  // projeta->SetSumw2();
  //TProfile* etaprofile = projeta->ProfileX("dNdeta_profile");
  //TProfile* etaprofile = projeta->ProfileX("dNdeta_profile");
  // TH1* etaprofile = projeta->ProjectionX("dNdeta_profile");
  //TProfile* etaprofile = fDataObject->ProfileX();
  /*
  TCanvas* ctest = new TCanvas("test","test",1200,600);
  ctest->Divide(3);
  ctest->cd(1);
  fDataObject->DrawCopy("box");
  ctest->cd(2);
  projeta->DrawCopy("colz");
  ctest->cd(3);
  etaprofile->DrawCopy();
  */
  //sumMult->Rebin(rebin);
  TH1F* hReldif = (TH1F*)sumMult->Clone("hReldif");
  if(rebin != 1) {
    RebinHistogram(sumMult,rebin);
    if(!realdata) {
      RebinHistogram(primMult,rebin);
      RebinHistogram(hReldif,rebin);
    }
  }
  hReldif->Add(primMult,-1);
  hReldif->Divide(primMult);
  TCanvas* c2 = new TCanvas("dN/deta","dN/deta",640,960);
  c2->Divide(1,2);//,0,0);
  c2->cd(2);
  gPad->Divide(1,2);
  gPad->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hReldif->SetTitle("");
  hReldif->GetYaxis()->SetTitle("Relative difference");
  hReldif->GetYaxis()->SetRangeUser(-0.2,0.2);
  hReldif->SetStats(0000);
  hReldif->GetXaxis()->SetTitle("#eta");
  
  hReldif->DrawCopy();
  //SPD Rel dif
  TH1F* hReldifSPD = (TH1F*)hSPDcombi->Clone("hReldifSPD");
  if(rebin != 1) {
    RebinHistogram(hSPDcombi,rebin);
    if(!realdata) {
      RebinHistogram(hReldifSPD,rebin);
    }
  }
  hReldifSPD->Add(primMult,-1);
  hReldifSPD->Divide(primMult);
  hReldifSPD->DrawCopy("same");
  
  TLine* zeroLine = new TLine(-4,0,6,0);
  zeroLine->SetLineWidth(2);
  zeroLine->SetLineColor(kBlack);
  zeroLine->Draw();
  hPrim->SetTitle("");
  hPrim->SetLabelFont(132,"xyz");
  hPrim->SetFillColor(kGray);
  primMult->SetFillColor(kBlue);
  
  
  c2->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hPrim->GetYaxis()->SetRangeUser(2,10);
  //hPrim->GetYaxis()->SetRangeUser(0,20);
  hPrim->SetStats(0000);
  hPrim->GetXaxis()->SetTitle("#eta");
  sumMult->SetTitle("");
  sumMult->SetStats(kFALSE);
  sumMult->SetXTitle("#eta");
  sumMult->SetYTitle("#frac{dN}{d#eta_{ch}}");
  // sumMult->GetYaxis()->SetRangeUser
  //sumMult->Scale(1/(Float_t)rebin);
  sumMult->DrawCopy("PE");
  //Syst errors
  TH1F* sumMultSystPos = (TH1F*)sumMult->Clone("systerrorsPos");
  TH1F* sumMultSystNeg = (TH1F*)sumMult->Clone("systerrorsNeg");
  for(Int_t jj = 1;jj <=sumMultSystPos->GetNbinsX();jj++) {
    if(sumMultSystPos->GetBinCenter(jj) < 0) {
      sumMultSystPos->SetBinError(jj,0);
      sumMultSystPos->SetBinContent(jj,0);
      continue;
    }
      
    if(sumMultSystPos->GetBinContent(jj) > 0)
      sumMultSystPos->SetBinError(jj,TMath::Sqrt(TMath::Power(sumMultSystPos->GetBinError(jj),2)+TMath::Power(0.1*sumMultSystPos->GetBinContent(jj),2)));
  }
  for(Int_t jj = 1;jj <=sumMultSystNeg->GetNbinsX();jj++) {
    if(sumMultSystNeg->GetBinCenter(jj) > 0) {
      sumMultSystNeg->SetBinError(jj,0);
      sumMultSystNeg->SetBinContent(jj,0);
      continue;
    }
      
    if(sumMultSystNeg->GetBinContent(jj) > 0)
      sumMultSystNeg->SetBinError(jj,TMath::Sqrt(TMath::Power(sumMultSystNeg->GetBinError(jj),2)+TMath::Power(0.1*sumMultSystNeg->GetBinContent(jj),2)));
  }
  
  
  
  sumMultSystPos->SetFillColor(18);
  sumMultSystNeg->SetFillColor(18);
  //sumMultSystPos->DrawCopy("sameE5");
  //sumMultSystNeg->DrawCopy("sameE5");
  //End syst errors
  
  
  
  if(what != kHits && what != kHitsTrVtx && hPrim->GetEntries())
    hPrim->DrawCopy("E3same");
  if(what == kHits || what == kHitsTrVtx)
    primMult->DrawCopy("E3same");
  hPrim->GetXaxis()->SetTitle("#eta");
  
  TH1F* hMirror = new TH1F("hMirror","mirrored",sumMult->GetNbinsX(),sumMult->GetXaxis()->GetXmin(),sumMult->GetXaxis()->GetXmax());
  
  for(Int_t i= (Int_t)0.5*sumMult->GetNbinsX(); i<=sumMult->GetNbinsX(); i++) {
    Float_t eta   = sumMult->GetBinCenter(i);
    Float_t value = sumMult->GetBinContent(i);
    Float_t error = sumMult->GetBinError(i);
    if(value > 0) {
      hMirror->SetBinContent(hMirror->FindBin(-1*eta),value);
      hMirror->SetBinError(hMirror->FindBin(-1*eta),error);
    }
  }
  TH1F* hReldifLR = (TH1F*)sumMult->Clone("hReldifLeftRight");
  hReldifLR->Add(hMirror,-1);
  hReldifLR->Divide(hMirror);
  c2->cd(2);
  gPad->cd(1);
  hReldifLR->GetYaxis()->SetRangeUser(-0.2,0.2);

  hReldifLR->SetYTitle("Relative difference left-right");
  hReldifLR->SetLabelSize(2*hReldifLR->GetLabelSize("X"),"X"); 
  hReldifLR->SetLabelSize(2*hReldifLR->GetLabelSize("Y"),"Y");
  hReldifLR->SetTitleSize(2*hReldifLR->GetTitleSize("X"),"X");
  hReldifLR->SetTitleSize(2*hReldifLR->GetTitleSize("Y"),"Y");
  hReldifLR->SetTitleOffset(0.7*hReldifLR->GetTitleOffset("X"),"X");
  hReldifLR->SetTitleOffset(0.5*hReldifLR->GetTitleOffset("Y"),"Y");
  
  
  if(realdata) 
    hReldifLR->DrawCopy();
  zeroLine->Draw();
  c2->cd(1);
  hMirror->SetMarkerStyle(25);
  
  hPrim->SetMarkerStyle(8);
  
  //graph2->Draw("P");
  TH1F* hPythiaMB = 0;
  
  TFile* fpyt = 0;
  
  if(!pythiafile.Contains("none"))
    fpyt = TFile::Open(pythiafile);
  
  if(realdata ) {
    if(fpyt && !pythiafile.Contains("none")) {
      hPythiaMB = (TH1F*)fpyt->Get("hPythiaMB");
    }
    else
      std::cout<<"no pythia for this energy"<<std::endl;
  }
  if(what != kHits && what != kHitsTrVtx) {
    if(hPrim->GetEntries() != 0 && !realdata)
      hPrim->DrawCopy("PE3same");
    if(realdata) {
      //graph->Draw("PEsame");
      //graph2->Draw("PEsame");
      if(what != kMultNSD) {
	graphinel->Draw("PEsame");
	graphinel2->Draw("PEsame");
	p7742D1x1y1->Draw("PEsame");
      }
      else{
	graph->Draw("PEsame");
	graph2->Draw("PEsame");
	p7742D2x1y1->Draw("PEsame");
	p7743D8x1y1->Draw("PEsame");

      }
	
    }

  }
  
  
  sumMult->DrawCopy("PEsame");
  
  
	
  hSPDcombi->DrawCopy("same");
  if(realdata)
    hMirror->DrawCopy("PEsame");
    
  TLegend* leg = new TLegend(0.3,0.2,0.7,0.45,"");
  if(!realdata) {
    if(what != kHits && what != kHitsTrVtx)
      leg->AddEntry(hPrim,"Primary","pf");
    else
      leg->AddEntry(primMult,"Hits","pf");
  }
  
  if(what == kMult)
    leg->AddEntry(sumMult,"FMD INEL","p");
  else if(what == kMultTrVtx)
    leg->AddEntry(sumMult,"FMD TrVtx","p");
  else if(what == kMultNSD)
    leg->AddEntry(sumMult,"FMD NSD","p");
  
  
  if(realdata) {
  if(what == kMult)
    leg->AddEntry(hMirror,"Mirror FMD INEL","p");
  else if(what == kMultTrVtx)
    leg->AddEntry(hMirror,"Mirror FMD TrVtx","p");
  else if(what == kMultNSD)
    leg->AddEntry(hMirror,"Mirror FMD NSD","p");
   
  if(what == kMultNSD) {  
    leg->AddEntry(graph,"UA5 NSD","p");
    leg->AddEntry(graph2,"Mirror UA5 NSD","p"); }
  else if(what == kMult) {
    leg->AddEntry(graphinel,"UA5 INEL","p");
    leg->AddEntry(graphinel2,"Mirror UA5 INEL","p");
  }
    //leg->AddEntry(hPythiaMB,"Pythia MB","l");
  leg->AddEntry(hSPDcombi,"HHD SPD clusters","p");
  if(what == kMult)
    leg->AddEntry(p7742D1x1y1,"ALICE INEL published","p");
  if(what == kMultNSD) {
    leg->AddEntry(p7742D2x1y1,"ALICE NSD published","p");
    leg->AddEntry(p7743D8x1y1, "CMS NSD published","p");
  }
  }
  leg->Draw();
  
  c2->cd(2);
  gPad->cd(2);
  TH1F* hRatioMultPythia = 0;
  TH1F* hRatioMultUA5 = 0;
  TH1F* hRatioMultUA5_rebin5 = new TH1F("hRatioMultUA5_rebin5","hRatioMultUA5_rebin5",sumMult->GetNbinsX(),sumMult->GetXaxis()->GetXmin(),sumMult->GetXaxis()->GetXmax());
  
  Double_t xval=0, yval=0;
  Int_t ipos=0;
  for(Int_t bb =hRatioMultUA5_rebin5->FindBin(0); bb<=hRatioMultUA5_rebin5->GetNbinsX();bb++) {
    
    
    
    xval=yval=0;
    
    if(hRatioMultUA5_rebin5->GetBinCenter(bb) >= 0) {
      if(what == kMultNSD)
	graph->GetPoint(ipos,xval,yval);
      else
	graphinel->GetPoint(ipos,xval,yval);
      
      if(yval>0) {
	hRatioMultUA5_rebin5->SetBinContent(bb,yval);
	hRatioMultUA5_rebin5->SetBinError(bb,graphinel->GetErrorY(ipos));
	if(hRatioMultUA5_rebin5->GetBinCenter(bb) < 4) {
	  hRatioMultUA5_rebin5->SetBinContent(hRatioMultUA5_rebin5->FindBin(-1*hRatioMultUA5_rebin5->GetBinCenter(bb)),yval);
	  hRatioMultUA5_rebin5->SetBinError(hRatioMultUA5_rebin5->FindBin(-1*hRatioMultUA5_rebin5->GetBinCenter(bb)),(what == kMultNSD ? graph->GetErrorY(ipos) : graph->GetErrorY(ipos)));
	
	}
	ipos++;
      }
    }
  }
       
      

  if(realdata) {
    
    Float_t ratio = 1;
    if(hPythiaMB) {
      if(hPythiaMB->GetNbinsX() !=  sumMult->GetNbinsX())
	ratio = (Float_t)sumMult->GetNbinsX() / (Float_t)hPythiaMB->GetNbinsX();
    }
    
    hRatioMultPythia = (TH1F*)sumMult->Clone("hRatioMultPythia");
    hRatioMultUA5    = (TH1F*)sumMult->Clone("hRatioMultUA5");
    if(ratio > 1) {
      hRatioMultPythia->Rebin((Int_t)ratio);
      hRatioMultPythia->Scale(1/ratio);
    }
    if(ratio < 1 && hPythiaMB) {
      hPythiaMB->Rebin((Int_t)(1/ratio));
      hPythiaMB->Scale(ratio);
    }
    
    if(rebin !=5) {
      TGraphErrors* tmp1;
      if(what == kMultNSD)
      tmp1 = (TGraphErrors*)graph->Clone("UA5tmp");
    else
      tmp1 = (TGraphErrors*)graphinel->Clone("UA5tmp");
    
    tmp1->Fit("pol8","Q0","Q0",1.5,6);
      TF1* hFit = tmp1->GetFunction("pol8");
      for(Int_t ii = hRatioMultUA5->GetNbinsX() / 2; ii<hRatioMultUA5->GetNbinsX();ii++) { 
	if(hRatioMultUA5->GetBinContent(ii) > 0) {
	  Float_t errorratio = hRatioMultUA5->GetBinError(ii) / hRatioMultUA5->GetBinContent(ii);
	  hRatioMultUA5->SetBinContent(ii,
				       hRatioMultUA5->GetBinContent(ii) / 
				       ((hRatioMultUA5->GetBinCenter(ii) > 4.8) ? hFit->Eval(hRatioMultUA5->GetBinCenter(ii-1)) : hFit->Eval(hRatioMultUA5->GetBinCenter(ii))));
	  hRatioMultUA5->SetBinError(ii,hRatioMultUA5->GetBinContent(ii) * TMath::Sqrt(0.02*0.02 + TMath::Power(errorratio,2)));
	  
	    
	}
	  
      }
      
      TGraphErrors* tmp2;
      if(what == kMultNSD)
	tmp2 = (TGraphErrors*)graph2->Clone("UA5tmp2");
      else 
	tmp2 = (TGraphErrors*)graphinel2->Clone("UA5tmp2");

      //tmp2->Fit("pol8","Q0","Q0",-3.7,-1.5);
      tmp2->Fit("pol7","Q0","Q0",-3.7,0);
      hFit = tmp2->GetFunction("pol7");
      for(Int_t ii = 0; ii<hRatioMultUA5->GetNbinsX()/2;ii++) { 
	if(hRatioMultUA5->GetBinContent(ii) > 0) {
	  Float_t errorratio = hRatioMultUA5->GetBinError(ii) / hRatioMultUA5->GetBinContent(ii);
	  hRatioMultUA5->SetBinContent(ii,hRatioMultUA5->GetBinContent(ii) / hFit->Eval(hRatioMultUA5->GetBinCenter(ii))  );
	  hRatioMultUA5->SetBinError(ii,hRatioMultUA5->GetBinContent(ii) * TMath::Sqrt(0.02*0.02 + TMath::Power(errorratio,2)));
	  
	}
	
	
	graphinel->GetHistogram()->SetStats(kFALSE);
	graphinel2->GetHistogram()->SetStats(kFALSE);
	
      }
    }
    else hRatioMultUA5->Divide(hRatioMultUA5_rebin5);
      //}
  
      
  
      gPad->SetGridx();
      gPad->SetGridy();
    TLine* oneLine = new TLine(-4,1,6,1);
    oneLine->SetLineWidth(2);
    
    hRatioMultUA5->SetMarkerColor(kGreen);
    hRatioMultUA5->SetMarkerStyle(20);
    hRatioMultUA5->GetYaxis()->SetRangeUser(0.75,1.25);
    
    hRatioMultUA5->SetYTitle("Ratio FMD/UA5");
    hRatioMultUA5->SetLabelSize(2*hRatioMultUA5->GetLabelSize("X"),"X"); 
    hRatioMultUA5->SetLabelSize(2*hRatioMultUA5->GetLabelSize("Y"),"Y");
    hRatioMultUA5->SetTitleSize(2*hRatioMultUA5->GetTitleSize("X"),"X");
    hRatioMultUA5->SetTitleSize(2*hRatioMultUA5->GetTitleSize("Y"),"Y");
    hRatioMultUA5->SetTitleOffset(0.7*hRatioMultUA5->GetTitleOffset("X"),"X");
    hRatioMultUA5->SetTitleOffset(0.5*hRatioMultUA5->GetTitleOffset("Y"),"Y");
    
    
    
    hRatioMultUA5->DrawCopy();
    //hRatioMultPythia->DrawCopy("same");
    oneLine->Draw("same");
  }
  
  TCanvas* cPaper = new TCanvas("FMD_SPD_UA5","FMD_SPD_UA5",800,600);
  cPaper->cd();
  
  sumMult->DrawCopy();
  //sumMultSystPos->DrawCopy("sameE5");
  // sumMultSystNeg->DrawCopy("sameE5");
  //hTestdNdeta->DrawCopy("same");
  //hSPDcombi->DrawCopy("same");
  if(what == kMult)
    p7742D1x1y1->Draw("PEsame");
  if(what == kMultNSD) {
    p7742D2x1y1->Draw("PEsame");
    p7743D8x1y1->Draw("PEsame");
  }

  TLegend* leg2 = new TLegend(0.3,0.2,0.7,0.45,"");
  if(what == kMult)
    leg2->AddEntry(sumMult,"FMD INEL","p");
  else if(what == kMultTrVtx)
    leg2->AddEntry(sumMult,"FMD TrVtx","p");
  else if(what == kMultNSD) 
    leg2->AddEntry(sumMult,"FMD NSD","p");
  
  if(realdata) {
    
    /* if(what == kMult)
    leg2->AddEntry(hMirror,"Mirror FMD INEL","p");
  else if(what == kMultTrVtx)
    leg2->AddEntry(hMirror,"FMD TrVtx","p");
  else if(what == kMultNSD)
  leg2->AddEntry(hMirror,"FMD NSD","p");*/
   
  if(what == kMultNSD) {  
    leg2->AddEntry(graph,"UA5 NSD","p");
    leg2->AddEntry(graph2,"Mirror UA5 NSD","p"); }
  else if(what == kMult) {
    leg2->AddEntry(graphinel,"UA5 INEL","p");
    leg2->AddEntry(graphinel2,"Mirror UA5 INEL","p");
  }
    //leg2->AddEntry(hPythiaMB,"Pythia MB","l");
    //leg2->AddEntry(hSPDcombi,"HHD SPD clusters","p");
  if(what == kMult)
    leg2->AddEntry(p7742D1x1y1,"ALICE INEL published","p");
  if(what == kMultNSD) {
    leg2->AddEntry(p7742D2x1y1,"ALICE NSD published","p");
    leg2->AddEntry(p7743D8x1y1,"CMS NSD published","p");
  }
  }
  leg2->Draw();
  
  
  
  if(what != kMultNSD) {
    graphinel->Draw("PEsame");
    graphinel2->Draw("PEsame");
  }
  else {
    graph->Draw("PEsame");
    graph2->Draw("PEsame");
  }
    
  sumMult->DrawCopy("PEsame");
  c2->cd(1);

  c2->Print("fmdana.png");
  TString species;
  
  switch(what) {
  case kMult: 
    species = "mult";
    break;
  case kMultTrVtx:
    species = "mult_TrVtx";
    break;
  case kHits:
    species = "hits";
    break;
  case kHitsTrVtx:
    species = "hits_TrVtx";
    break;
  case kMultNSD:
    species = "mult_NSD";
    break;
    
  default:
    AliWarning("no valid Analysis entry!");
    break;
  }
  TFile fout(Form("fmd_dNdeta_%s.root",species.Data()),"recreate");
  if(hRatioMultPythia)
    hRatioMultPythia->Write();
  hPrim->Write();
  graph->Write("UA5nsd");
  graph2->Write("UA5mirrornsd");
  graphinel->Write("UA5inel");
  graphinel2->Write("UA5mirrorinel");
  sumMult->Write();
  hSPDcombi->Write();
  hReldif->Write();
  fMultList[what]->Write();
  c2->Write();
  //fDataObject->Write();
  fout.Close();
    
}
//_____________________________________________________________________
void AliFMDDndeta::CreateSharingEfficiency(const Char_t* filename, Bool_t store) {
  //Generate sharing efficiency
  Init(filename);
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  //pars->Init(kTRUE,AliFMDAnaParameters::kBackgroundCorrection);
  Int_t nVertexBins = pars->GetNvtxBins();
  
  SetNames(kHits);
  // "nEvents";
  TH1I* hEvents          = (TH1I*)fList->FindObject(fEvents.Data());
  // "nMCEventsNoCuts"
  TH1I* hPrimEvents      = (TH1I*)fList->FindObject(fPrimEvents.Data());

  SetNames(kHitsTrVtx);
  // "nEvents";
  TH1I* hEventsTrVtx     = (TH1I*)fList->FindObject(fEvents.Data());
  // "nMCEvents"
  TH1I* hPrimEventsTrVtx = (TH1I*)fList->FindObject(fPrimEvents.Data());
  
  AliFMDAnaCalibSharingEfficiency* sharEff = new AliFMDAnaCalibSharingEfficiency();
  
  for(Int_t det = 1; det<=3; det++) {
    Int_t maxRing = (det == 1 ? 0 : 1);
    for(Int_t ring = 0;ring<=maxRing;ring++) {
      Char_t ringChar = (ring == 0 ? 'I' : 'O');
      for(Int_t v=0; v< nVertexBins; v++) {
	// Get histograms like  hits_NoCuts_FMD%d%c_vtxbin%d_proj
	// - from AliFMDAnalysisTaskBackgroundCorrection.cxx
	TH1F* hHits = (TH1F*)fList->FindObject(GetAnalysisName(kHits,det, ringChar, v));
	// Get histograms like hits_FMD%d%c_vtxbin%d_proj 
	// - from AliFMDAnalysisTaskBackgroundCorrection.cxx 
	TH1F* hHitsTrVtx   = (TH1F*)fList->FindObject(GetAnalysisName(kHitsTrVtx,det, ringChar, v));
	// Get histograms like "hMCHits_nocuts_FMD%d%c_vtxbin%d"
	// - from AliFMDAnalysisTaskSharing.cxx
	TH1F* hMCHits      = (TH1F*)fList->FindObject(GetPrimName(kHits,det, ringChar, v));
	// Get histograms like "hMCHits_FMD%d%c_vtxbin%d"
	// - from AliFMDAnalysisTaskSharing.cxx
	TH1F* hMCHitsTrVtx = (TH1F*)fList->FindObject(GetPrimName(kHitsTrVtx,det, ringChar, v));
	
	TH1F* hCorrection  = (TH1F*)hHits->Clone(Form("hCorrection_FMD%d%c_vtx%d",det, ringChar, v));
	TH1F* hCorrectionTrVtx  = (TH1F*)hHitsTrVtx->Clone(Form("hCorrection_FMD%d%c_vtx%d",det, ringChar, v));
	if(hEvents->GetBinContent(v+1))
	  hCorrection->Scale(1./(Float_t)hEvents->GetBinContent(v+1));
	if(hPrimEvents->GetBinContent(v+1))
	  hMCHits->Scale(1./(Float_t)hPrimEvents->GetBinContent(v+1));
	if(hEventsTrVtx->GetBinContent(v+1))
	  hCorrectionTrVtx->Scale(1./(Float_t)hEventsTrVtx->GetBinContent(v+1));
	if(hPrimEventsTrVtx->GetBinContent(v+1))
	  hMCHitsTrVtx->Scale(1./(Float_t)hPrimEventsTrVtx->GetBinContent(v+1));
	
	hCorrection->Divide(hMCHits);
	hCorrectionTrVtx->Divide(hMCHitsTrVtx);
	
	//sharEff->SetSharingEff(det,ringChar,v,hCorrection);
	sharEff->SetSharingEff(det,ringChar,v,hCorrection);
	sharEff->SetSharingEffTrVtx(det,ringChar,v,hCorrectionTrVtx);
	
	//	std::cout<<hHits<<"  "<<hHitsTrVtx<<"   "<<hPrim<<"    "<<hPrimTrVtx<<std::endl;

	}
      }
    }

  if(store) {
    TFile fsharing(pars->GetPath(pars->GetSharingEffID()),"RECREATE");
    sharEff->Write(AliFMDAnaParameters::GetSharingEffID());
    fsharing.Close();
  }
  
}
//_____________________________________________________________________
void AliFMDDndeta::RebinHistogram(TH1F* hist, Int_t rebin) {
  //Clever rebinner
  
  if(hist->GetNbinsX()%rebin) {
    std::cout<<"cannot rebin "<<hist->GetNbinsX()%rebin<< "is not zero"<<std::endl;
    return;

  }
  TH1F* cloneHist = (TH1F*)hist->Clone();
  cloneHist->Rebin(rebin);
  
  Int_t nBinsNew = hist->GetNbinsX() / rebin;
  
  for(Int_t i =1;i<=nBinsNew; i++) {
    Float_t bincontent = 0;
    Float_t sumofw = 0;
    Float_t sumformean = 0;
    Int_t   nBins = 0;
    for(Int_t j=1; j<=rebin;j++) {
      if(hist->GetBinContent((i-1)*rebin + j) > 0) {
	bincontent = bincontent + hist->GetBinContent((i-1)*rebin + j);
	nBins++;
	Float_t weight = 1 / TMath::Power(hist->GetBinError((i-1)*rebin + j),2);
	sumofw = sumofw + weight;
	sumformean = sumformean + weight*hist->GetBinContent((i-1)*rebin + j);
	
      }
    }
    
    if(bincontent > 0 ) {
      cloneHist->SetBinContent(i,sumformean / sumofw);
      cloneHist->SetBinError(i,TMath::Sqrt(sumformean));//TMath::Sqrt(1/sumofw));
    }
  }
  hist->Rebin(rebin);
  for(Int_t i =1;i<=nBinsNew; i++) {
    hist->SetBinContent(i,cloneHist->GetBinContent(i));
  }
  
  
}
//_____________________________________________________________________
//
// EOF
//
// EOF
