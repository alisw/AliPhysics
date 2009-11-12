
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
#include "TPad.h"
#include "iostream"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
ClassImp(AliFMDDndeta)
//_____________________________________________________________________

AliFMDDndeta::AliFMDDndeta() 
: TObject(),
  fList(0),
  fMultList(),
  fNbinsToCut(2),
  fVtxCut(10),
  fIsInit(kFALSE),
  fIsGenerated(),
  fPrimEvents(),
  fEvents(),
  fPrimdNdeta()
{
 
}
//_____________________________________________________________________
void AliFMDDndeta::SetNames(Analysis what) {
  
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
  default:
    break;
  }
}
//_____________________________________________________________________
const char* AliFMDDndeta::GetAnalysisName(Analysis what, UShort_t det, Char_t ring, Int_t vtxbin) {
  
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
  default:
    break;
  }
  //std::cout<<name.Data()<<std::endl;
  return name;
}
//_____________________________________________________________________
const char* AliFMDDndeta::GetPrimName(Analysis what, UShort_t det, Char_t ring, Int_t vtxbin) {
  
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
  default:
    break;
  }
  return name;
}
//_____________________________________________________________________
void AliFMDDndeta::GenerateHits() {

  AliFMDAnaParameters* pars =  AliFMDAnaParameters::Instance();
  TH2F* hTmp         =  pars->GetBackgroundCorrection(1,'I',0);
  Int_t nVertexBins  =  pars->GetNvtxBins();
  Int_t nEtaBins     =  hTmp->GetNbinsX();
  Float_t etaMin     =  hTmp->GetXaxis()->GetXmin(); 
  Float_t etaMax     =  hTmp->GetXaxis()->GetXmax();

  for(Int_t i = 0; i<nVertexBins;i++) {
    TH1F* hHits = new TH1F(Form("hMCHits_vtxbin%d",i),Form("hHits_vtxbin%d",i),nEtaBins,etaMin,etaMax);
    hHits->Sumw2();
    fMultList.Add(hHits);
  }
   
  for(Int_t det = 1; det<=3; det++) {
    Int_t maxRing = (det == 1 ? 0 : 1);
    for(Int_t ring = 0;ring<=maxRing;ring++) {
      
      Char_t ringChar = (ring == 0 ? 'I' : 'O');
      for(Int_t v=0; v< nVertexBins; v++) {
	TH1F* hits = (TH1F*)fList->FindObject(GetPrimName(kHits,det,ringChar,v));
	
	Int_t nNonZero = 0, nNonZeroInData = 0;
	
	//removing edges
	for(Int_t i =1;i<=hits->GetNbinsX();i++) {
	  if(hits->GetBinContent(i) !=0)
	    nNonZero++;
	}

	TH1F* sumMultHist = (TH1F*)fMultList.FindObject(Form("hMCHits_vtxbin%d",v));
	for(Int_t i =1;i<=sumMultHist->GetNbinsX();i++) {
	  if(hits->GetBinContent(i) == 0 ) continue;
	  
	  nNonZeroInData++;
	  
	  if(nNonZeroInData <=fNbinsToCut || nNonZeroInData > (nNonZero - fNbinsToCut)) {
	    continue;
	  }
	  if(sumMultHist->GetBinContent(i) == 0 && hits->GetBinContent(i) != 0){
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
  
  TFile* fin = TFile::Open(filename);
  
  if(!fin) {
    AliWarning("No file - aborting !");
    return;
  }
  AliFMDAnaParameters* pars =  AliFMDAnaParameters::Instance();
  //pars->Init();
  
  fList = (TList*)fin->Get(Form("%s/BackgroundCorrected",pars->GetDndetaAnalysisName()));
  
  if(!fList) //an old file ? Perhaps...
    fList = (TList*)fin->Get("BackgroundCorrected");
    
  fIsGenerated[kHits]      = kFALSE;
  fIsGenerated[kMult]      = kFALSE;  
  fIsGenerated[kMultTrVtx] = kFALSE;
  
  fIsInit = kTRUE;
}
//_____________________________________________________________________
void AliFMDDndeta::GenerateMult(Analysis what) {
 
  
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
  
  if(what == kHits)
    GenerateHits();
  
  TH1I* hMCEvents = (TH1I*)fList->FindObject(fPrimEvents.Data());
  
  AliFMDAnaParameters* pars =  AliFMDAnaParameters::Instance();
  
  TH2F* hTmp         =  pars->GetBackgroundCorrection(1,'I',0);
  Int_t nVertexBins  =  pars->GetNvtxBins();
  Int_t nEtaBins     =  hTmp->GetNbinsX();
  Float_t etaMin     =  hTmp->GetXaxis()->GetXmin(); 
  Float_t etaMax     =  hTmp->GetXaxis()->GetXmax();
  
  for(Int_t i = 0; i<nVertexBins;i++) {
    TH1F* hMult = new TH1F(Form("hMult_vtxbin%d",i),Form("hMult_vtxbin%d",i),nEtaBins,etaMin,etaMax);
    hMult->Sumw2();
    fMultList.Add(hMult);
  }
  
  for(Int_t det = 1; det<=3;det++) {
    Int_t maxRing = (det == 1 ? 0 : 1);
    for(Int_t iring = 0; iring<=maxRing; iring++) {
      Char_t ringChar = (iring == 0 ? 'I' : 'O');
      TH1F* hRingMult= new TH1F(Form("hRingMult_FMD%d%c",det,ringChar),Form("hRingMult_FMD%d%c",det,ringChar),nEtaBins,etaMin,etaMax);
      fMultList.Add(hRingMult);
    }
  }
  TH1I* hEvents = (TH1I*)fList->FindObject(fEvents.Data());
  TH1F*  hPrimVtx = 0;
  
  	
    
  for(Int_t det = 1; det<=3; det++) {
    Int_t maxRing = (det == 1 ? 0 : 1);
    for(Int_t ring = 0;ring<=maxRing;ring++) {
      
      Char_t ringChar = (ring == 0 ? 'I' : 'O');
      for(Int_t v=0; v< nVertexBins; v++) {
	if(det == 1) {
	  if(what == kHits)
	    hPrimVtx = (TH1F*)fMultList.FindObject(Form("hMCHits_vtxbin%d",v));
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
	
	if(vtxZ1<(-1*fVtxCut) || vtxZ2 >fVtxCut)
	  continue;
	Float_t nEvents = (Float_t)hEvents->GetBinContent(v+1);
	
	//TH1F* multhistproj = (TH1F*)fList->FindObject(Form("dNdeta_FMD%d%c_vtxbin%d_proj",det,ringChar,v));
       	TH1F* multhistproj = (TH1F*)fList->FindObject(GetAnalysisName(what,det,ringChar,v));
	if(nEvents)
	  multhistproj->Scale(1/nEvents);
	Float_t xnBins = multhistproj->GetNbinsX();
	Float_t xrange = multhistproj->GetXaxis()->GetXmax() - multhistproj->GetXaxis()->GetXmin(); 
	multhistproj->Scale(xnBins / xrange);
	

	
	Int_t nNonZero = 0, nNonZeroInData = 0;
	
	//removing edges
	
	for(Int_t i =1;i<=multhistproj->GetNbinsX();i++) {
	  if(multhistproj->GetBinContent(i) !=0)
	    nNonZero++;
	}
	Int_t nBinsOld = fNbinsToCut;
	//	if(det == 2 && ringChar =='I') {
	//  fNbinsToCut = 1;
	//	}
	TH1F* hRingMult = (TH1F*)fMultList.FindObject(Form("hRingMult_FMD%d%c",det,ringChar));
	
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
	
	TH1F* sumMultHist = (TH1F*)fMultList.FindObject(Form("hMult_vtxbin%d",v));
	
	for(Int_t i =1;i<=sumMultHist->GetNbinsX();i++) {
	  if(multhistproj->GetBinContent(i) != 0 ) {	  
	    nNonZeroInData++;
	  
	    if(nNonZeroInData <=fNbinsToCut || nNonZeroInData > (nNonZero - fNbinsToCut)) {
	      continue;
	    }
	    if(sumMultHist->GetBinContent(i) == 0 && multhistproj->GetBinContent(i) != 0){
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
  
  TH1F* sumMult  = new TH1F("hSumMult","hSumMult",nEtaBins,etaMin,etaMax);
  sumMult->Sumw2();
  fMultList.Add(sumMult);
  TH1F* primMult  = new TH1F("hPrimMult","hPrimMult",nEtaBins,etaMin,etaMax);
  primMult->Sumw2();
  fMultList.Add(primMult);

  Float_t wav  = 0, sumofw=0, weight = 0, primwav  = 0, primsumofw=0, primweight = 0;//, etaofbin = 0;
  for(Int_t i =1; i<=sumMult->GetNbinsX();i++) {
    
    for(Int_t v = 0; v<nVertexBins;v++) {
      if(what == kHits)
	hPrimVtx = (TH1F*)fMultList.FindObject(Form("hMCHits_vtxbin%d",v));
      else
	hPrimVtx = (TH1F*)fList->FindObject(GetPrimName(what,0,0,v));
      TH1F* sumMultHist = (TH1F*)fMultList.FindObject(Form("hMult_vtxbin%d",v));
      //etaofbin += sumMultHist->GetBinCenter(i);
      // if( hPrimVtx->GetBinContent(i)!=0 && hPrimVtx->GetBinError(i)>0.001 && sumMultHist->GetBinContent(i)!=0) {
      if( hPrimVtx->GetBinContent(i)!=0) {
	primweight = 1/TMath::Power(hPrimVtx->GetBinError(i),2);
	primwav    = primwav + primweight*hPrimVtx->GetBinContent(i);
	primsumofw = primsumofw + primweight;
      }
      if(sumMultHist->GetBinContent(i)!=0) {
	
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
void AliFMDDndeta::DrawDndeta(Analysis what, Int_t rebin) {
  
  if(!fIsInit) {
    AliWarning("Not initialised - call Init to remedy");
    return;
  }

  if(!fIsGenerated[what]) {
    AliWarning("Not generated - generate first then draw!");
    return;
  }
    
  SetNames(what);
  
  AliFMDAnaParameters* pars =  AliFMDAnaParameters::Instance();
  TCanvas* c1 = new TCanvas("cvertexbins","Overlaps",1400,800);
  c1->Divide(5,2);
  Int_t nVertexBins  =  pars->GetNvtxBins();
  
  TH1I* hEvents   = (TH1I*)fList->FindObject(fEvents.Data());
  TH1I* hMCEvents = (TH1I*)fList->FindObject(fPrimEvents.Data());
  
  TH1F*  hPrimVtx = 0;
  for(Int_t det = 1; det<=3; det++) {
    Int_t maxRing = (det == 1 ? 0 : 1);
    for(Int_t ring = 0;ring<=maxRing;ring++) {
      for(Int_t v=0; v< nVertexBins; v++) {
	Char_t ringChar = (ring == 0 ? 'I' : 'O');
	Double_t delta = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
	Float_t vtxZ1   = (delta*v) - pars->GetVtxCutZ();
	Float_t vtxZ2   = (delta*(v+1)) - pars->GetVtxCutZ();
	  
	if(vtxZ1<(-1*fVtxCut) || vtxZ2 >fVtxCut)
	  continue;
	Float_t nEvents = (Float_t)hEvents->GetBinContent(v+1);
		
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
	  if(what == kHits)
	    hPrimVtx = (TH1F*)fMultList.FindObject(Form("hMCHits_vtxbin%d",v));
	  else
	    hPrimVtx = (TH1F*)fList->FindObject(GetPrimName(what,det,ringChar,v));
	  
	  hPrimVtx->SetTitle("");
	  TLatex* l = new TLatex(0.14,0.92,Form("Vtx range [%.1f, %.1f], %d events",vtxZ1,vtxZ2,(Int_t)nEvents));
	  l->SetNDC(kTRUE);
	  hPrimVtx->SetFillColor(kGray);
	  hPrimVtx->SetLabelFont(132,"xyz");
	  hPrimVtx->SetStats(0000);
	  hPrimVtx->GetXaxis()->SetTitle("#eta");
	  hPrimVtx->GetYaxis()->SetRangeUser(0,15);
	  hPrimVtx->DrawCopy("E3");
	  l->Draw();
	  multhistproj->DrawCopy("same");
	  
	}
	else
	  multhistproj->DrawCopy("same");
	
	
      }
    }
  }
  
  for(Int_t v=0; v< nVertexBins; v++) {
    TH1F* sumMultHist = (TH1F*)fMultList.FindObject(Form("hMult_vtxbin%d",v));
    c1->cd(v+1);
    sumMultHist->SetMarkerStyle(25);
    sumMultHist->DrawCopy("same");
    
  }
  TH1F* primMult = (TH1F*)fMultList.FindObject("hPrimMult");
  TH1F* sumMult  = (TH1F*)fMultList.FindObject("hSumMult");
  sumMult->SetMarkerStyle(23);
  sumMult->SetMarkerColor(kRed);
  for(Int_t det = 1; det<=3;det++) {
    Int_t maxRing = (det == 1 ? 0 : 1);
    for(Int_t iring = 0; iring<=maxRing; iring++) {
      Char_t ringChar = (iring == 0 ? 'I' : 'O');
      TH1F* hRingMult= (TH1F*)fMultList.FindObject(Form("hRingMult_FMD%d%c",det,ringChar));
      hRingMult->Add(primMult,-1);
      hRingMult->Divide(primMult);
    }
  }
  
  
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
  
  /*
  TFile testin("/home/canute/ALICE/FMDanalysis/GridAnalysis/fmd_dNdeta_hits.root","READ");
    // TFile testin("/home/canute/ALICE/FMDanalysis/productionData/fmd_dNdeta_hits.root","READ");
  TH1F* hcorrect = (TH1F*)testin.Get("hReldif");
  hcorrect->SetName("djarnis");
  std::cout<<hcorrect<<std::endl;
  for(Int_t bb = 1;bb<=hcorrect->GetNbinsX();bb++) {
    hcorrect->SetBinContent(bb,hcorrect->GetBinContent(bb)+1);
  }
 
  sumMult->Divide(hcorrect);
  //delete hcorrect;
  */
  
  sumMult->Rebin(rebin);
  sumMult->Scale(1/rebin);
  primMult->Rebin(rebin);
  primMult->Scale(1/rebin);
  hPrim->Rebin(rebin);
  hPrim->Scale(1/rebin);
  if(what != kHits)
    primMult = hPrim;
  TH1F* hReldif = (TH1F*)sumMult->Clone("hReldif");
  hReldif->Add(primMult,-1);
  hReldif->Divide(primMult);
  //  hReldif->Add(hPrim,-1);
  // hReldif->Divide(hPrim);
  
  TCanvas* c2 = new TCanvas("dN/deta","dN/deta",640,960);
  c2->Divide(1,2);//,0,0);
  c2->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hReldif->SetTitle("");
  hReldif->GetYaxis()->SetTitle("Relative difference");
  hReldif->GetYaxis()->SetRangeUser(-0.2,0.2);
  hReldif->SetStats(0000);
  hReldif->GetXaxis()->SetTitle("#eta");
  hReldif->DrawCopy();
  TLine* zeroLine = new TLine(-6,0,6,0);
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
  sumMult->DrawCopy("PE");
  if(what != kHits)
    hPrim->DrawCopy("E3same");
  if(what == kHits)
    primMult->DrawCopy("E3same");
  hPrim->GetXaxis()->SetTitle("#eta");
  
  hPrim->SetMarkerStyle(8);
  if(what != kHits)
    hPrim->DrawCopy("PE3same");
  sumMult->DrawCopy("PEsame");
  std::cout<<"FMD 1 "<<sumMult->Integral(sumMult->FindBin(3),sumMult->FindBin(5),"width")<<std::endl;
  /*TFile* itsfile = TFile::Open(itsfilename);
  if(itsfile) {
    //itsfile->cd();
    
    //TList* itslist = (TList*)itsfile->Get("ITS_analysis");
    //TH1F* itshist0 = (TH1F*)itslist->FindObject("dndeta_check_0");
    //TH1F* itshist1 = (TH1F*)itslist->FindObject("dndeta_check_1");
    //TH1F* itshist2 = (TH1F*)itslist->FindObject("dndeta_check_2");
    TH1F* itshist0 = (TH1F*)itsfile->Get("dndeta/dNdEta");
    TH1F* itshist1= (TH1F*)itsfile->Get("dndeta/dNdEta_1");
    TH1F* itshist2 = (TH1F*)itsfile->Get("dndeta/dNdEta_2");
    itshist0->SetMarkerStyle(27);
    itshist1->SetMarkerStyle(28);
    itshist2->SetMarkerStyle(30);
    itshist0->DrawCopy("same");
    itshist1->DrawCopy("same");
    itshist2->DrawCopy("same");
  
  }
  */
  
  
  TLegend* leg = new TLegend(0.35,0.2,0.65,0.35,"");
  if(what != kHits)
    leg->AddEntry(hPrim,"Primary","pf");
  else
    leg->AddEntry(primMult,"Hits","pf");
  //leg->AddEntry(primMult,"Analysed prim","f");
  leg->AddEntry(sumMult,"Analysis","p");
  leg->Draw();
  TString species;
  
  switch(what) {
  case kMult: 
    species.Form("mult");
    break;
  case kMultTrVtx:
    species.Form("mult_TrVtx");
    break;
  case kHits:
    species.Form("hits");
    break;
  default:
    AliWarning("no valid Analysis entry!");
    break;
  }
  TFile fout(Form("fmd_dNdeta_%s.root",species.Data()),"recreate");
  sumMult->Write();
  hReldif->Write();
  fMultList.Write();
  fout.Close();
    
}
//_____________________________________________________________________
void AliFMDDndeta::CreateSharingEfficiency(const Char_t* filename, Bool_t store) {

  Init(filename);
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  //pars->Init(kTRUE,AliFMDAnaParameters::kBackgroundCorrection);
  Int_t nVertexBins = pars->GetNvtxBins();
  
  SetNames(kHits);
  TH1I* hEvents          = (TH1I*)fList->FindObject(fEvents.Data());
  TH1I* hPrimEvents      = (TH1I*)fList->FindObject(fPrimEvents.Data());

  SetNames(kHitsTrVtx);
  TH1I* hEventsTrVtx     = (TH1I*)fList->FindObject(fEvents.Data());
  TH1I* hPrimEventsTrVtx = (TH1I*)fList->FindObject(fPrimEvents.Data());
  
  AliFMDAnaCalibSharingEfficiency* sharEff = new AliFMDAnaCalibSharingEfficiency();
  
  for(Int_t det = 1; det<=3; det++) {
    Int_t maxRing = (det == 1 ? 0 : 1);
    for(Int_t ring = 0;ring<=maxRing;ring++) {
      Char_t ringChar = (ring == 0 ? 'I' : 'O');
      for(Int_t v=0; v< nVertexBins; v++) {
	TH1F* hHits = (TH1F*)fList->FindObject(GetAnalysisName(kHits,det, ringChar, v));
	TH1F* hHitsTrVtx   = (TH1F*)fList->FindObject(GetAnalysisName(kHitsTrVtx,det, ringChar, v));
	TH1F* hMCHits      = (TH1F*)fList->FindObject(GetPrimName(kHits,det, ringChar, v));
	TH1F* hMCHitsTrVtx = (TH1F*)fList->FindObject(GetPrimName(kHitsTrVtx,det, ringChar, v));
	
	TH1F* hCorrection  = (TH1F*)hHits->Clone(Form("hCorrection_FMD%d%c_vtx%d"));
	TH1F* hCorrectionTrVtx  = (TH1F*)hHitsTrVtx->Clone(Form("hCorrection_FMD%d%c_vtx%d"));
	
	hCorrection->Scale(1./(Float_t)hEvents->GetBinContent(v+1));
	hMCHits->Scale(1./(Float_t)hPrimEvents->GetBinContent(v+1));
	hCorrectionTrVtx->Scale(1./(Float_t)hEventsTrVtx->GetBinContent(v+1));
	hMCHitsTrVtx->Scale(1./(Float_t)hPrimEventsTrVtx->GetBinContent(v+1));
	hCorrection->Divide(hMCHits);
	hCorrectionTrVtx->Divide(hMCHitsTrVtx);
	
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
//
// EOF
//
