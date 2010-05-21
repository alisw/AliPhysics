 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>
#include <TMath.h>
//#include "AliFMDDebug.h"
#include "AliFMDAnalysisTaskSharing.h"
#include "AliAnalysisManager.h"
#include "AliESDFMD.h"
//#include "AliFMDGeometry.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliFMDAnaParameters.h"
#include "TH1F.h"
#include "TObjString.h"
//#include "/home/canute/ALICE/AliRoot/PWG0/AliPWG0Helper.h"
//#include "AliFMDParameters.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliMCParticle.h"
#include "AliFMDStripIndex.h"
#include "AliESDVZERO.h"
#include "AliESDtrack.h"

// This is the task to do the FMD sharing or hit merging.
// It reads the input ESDFMD data and posts an ESDFMD object to
// the tasks that must be performed after this task ie.
// Density, BackgroundCorrection and Dndeta.
// Author: Hans Hjersing Dalsgaard, hans.dalsgaard@cern.ch
 

ClassImp(AliFMDAnalysisTaskSharing)

//_____________________________________________________________________
AliFMDAnalysisTaskSharing::AliFMDAnalysisTaskSharing()
: fDebug(0),
  fESD(0x0),
  foutputESDFMD(),
  fSharedThis(kFALSE),
  fSharedPrev(kFALSE),
  fDiagList(0),
  fStandalone(kTRUE),
  fEsdVertex(0),
  fStatus(kTRUE),
  fLastTrackByStrip(0),
  fLastOrbit(0)
{
  // Default constructor
  DefineInput (0, AliESDEvent::Class());
  DefineOutput(0, AliESDFMD::Class());
  DefineOutput(1, AliESDVertex::Class());
  DefineOutput(2, AliESDEvent::Class());
  DefineOutput(3, TList::Class());
}
//_____________________________________________________________________
AliFMDAnalysisTaskSharing::AliFMDAnalysisTaskSharing(const char* name, Bool_t SE):
    AliAnalysisTask(name, "AnalysisTaskFMD"),
    fDebug(0),
    fESD(0x0),
    foutputESDFMD(),
    fSharedThis(kFALSE),
    fSharedPrev(kFALSE),
    fDiagList(0),
    fStandalone(kTRUE),
    fEsdVertex(0),
    fStatus(kTRUE),
    fLastTrackByStrip(0),
    fLastOrbit(0)
{
  // named constructor
  fStandalone = SE;
  if(fStandalone) {
    DefineInput (0, AliESDEvent::Class());
    DefineOutput(0, AliESDFMD::Class());
    DefineOutput(1, AliESDVertex::Class());
    DefineOutput(2, AliESDEvent::Class());
    DefineOutput(3, TList::Class());
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::CreateOutputObjects()
{
  // Create the output objects
  if(!foutputESDFMD)
    foutputESDFMD = new AliESDFMD();
  
  if(!fEsdVertex)
    fEsdVertex    = new AliESDVertex();
  //Diagnostics
  if(!fDiagList)
    fDiagList = new TList();
  
  fDiagList->SetName("Sharing diagnostics");
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  TH2F* hBg   = pars->GetBackgroundCorrection(1, 'I', 0);
  TH1F* hPrimary = new TH1F("hMultvsEtaNoCuts","hMultvsEtaNoCuts",
			    hBg->GetNbinsX(),
			    hBg->GetXaxis()->GetXmin(),
			    hBg->GetXaxis()->GetXmax());
  hPrimary->Sumw2();
  fDiagList->Add(hPrimary);
  TH1F* hXvtx = new TH1F("hXvtx","x vertex distribution",100,-2,2);
  TH1F* hYvtx = new TH1F("hYvtx","y vertex distribution",100,-2,2);
  TH1F* hZvtx = new TH1F("hZvtx","z vertex distribution",4*pars->GetNvtxBins(),-4*pars->GetVtxCutZ(),4*pars->GetVtxCutZ());
  
  fDiagList->Add(hXvtx);
  fDiagList->Add(hYvtx);
  fDiagList->Add(hZvtx);
  
  TH1F* hPrimVertexBin = 0;
  TH1F* hHits = 0;
  for(Int_t i = 0; i< pars->GetNvtxBins(); i++) {
    
    hPrimVertexBin = new TH1F(Form("primmult_NoCuts_vtxbin%d",i),
			      Form("primmult_NoCuts_vtxbin%d",i),
			      hBg->GetNbinsX(),
			      hBg->GetXaxis()->GetXmin(),
			      hBg->GetXaxis()->GetXmax());
    hPrimVertexBin->Sumw2();
    fDiagList->Add(hPrimVertexBin);
    
  }
  
  for(Int_t det = 1; det<=3; det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    
    for(Int_t iring = 0;iring<nRings; iring++) {
      Char_t ringChar = (iring == 0 ? 'I' : 'O');
      TH1F* hEdist        = new TH1F(Form("Edist_before_sharing_FMD%d%c", det, ringChar),
				     Form("Edist_before_sharing_FMD%d%c", det, ringChar),
				     1000,0,25);
      TH1F* hEdistAfter  = new TH1F(Form("Edist_after_sharing_FMD%d%c", det, ringChar),
				     Form("Edist_after_sharing_FMD%d%c", det, ringChar),
				     1000,0,25);
      
      
      //TH1F* hNstripsHit    = new TH1F(Form("N_strips_hit_FMD%d%c",det,ringChar),
      //				     Form("N_strips_hit_FMD%d%c",det,ringChar),
      //				     25,0,25);
      fDiagList->Add(hEdist);
      fDiagList->Add(hEdistAfter);
      //fDiagList->Add(hNstripsHit);
      
      for(Int_t i = 0; i< pars->GetNvtxBins(); i++) {
	hHits  = new TH1F(Form("hMCHits_nocuts_FMD%d%c_vtxbin%d",det,ringChar,i),Form("hMCHits_FMD%d%c_vtxbin%d",det,ringChar,i),
			  hBg->GetNbinsX(),
			  hBg->GetXaxis()->GetXmin(),
			  hBg->GetXaxis()->GetXmax());
	hHits->Sumw2();
	fDiagList->Add(hHits);

      }
      
    }
  }
  TH2F*  hCorrelationFMDSPDhits = new TH2F("hCorrelationFMDSPDhits","hCorrelationFMDSPDhits;SPD;FMD ",100,0,200,100,0,500);
  TH2F*  hCorrelationFMDSPD = new TH2F("hCorrelationFMDSPD","hCorrelationFMDSPD;SPD;FMD ",100,0,200,100,0,500);
  TH2F*  hCorrelationFMD1SPD = new TH2F("hCorrelationFMD1SPD","hCorrelationFMD1SPD;SPD;FMD1 ",100,0,200,100,0,200);
  TH2F*  hCorrelationFMD2ISPD = new TH2F("hCorrelationFMD2ISPD","hCorrelationFMD2ISPD;SPD;FMD2I ",100,0,200,100,0,200);
  TH2F*  hCorrelationFMD2OSPD = new TH2F("hCorrelationFMD2OSPD","hCorrelationFMD2OSPD;SPD;FMD2O ",100,0,200,100,0,200);
  TH2F*  hCorrelationFMD3ISPD = new TH2F("hCorrelationFMD3ISPD","hCorrelationFMD3ISPD;SPD;FMD3I ",100,0,200,100,0,200);
  TH2F*  hCorrelationFMD3OSPD = new TH2F("hCorrelationFMD3OSPD","hCorrelationFMD3OSPD;SPD;FMD3O ",100,0,200,100,0,200);
  TH2F*  hCorrelationFMDGoodtracks = new TH2F("hCorrelationFMDGoodtracks","hCorrelationGoodtracks;good tracks;FMD ",100,0,200,100,0,200);
  TH2F*  hCorrelationFMDBadtracks = new TH2F("hCorrelationFMDBadtracks","hCorrelationBadtracks;bad tracks;FMD ",100,0,200,100,0,200);
  TH2F*  hCorrelationGoodbadtracks = new TH2F("hCorrelationGoodbadtracks","hCorrelationGoodbadtracks;good tracks;bad tracks ",100,0,200,100,0,200);
  TH2F*  hCorrelationSPDTracklets = new TH2F("hCorrelationSPDTracklets","hCorrelationSPDTracklets;hits ; tracklets ",100,0,500,100,0,200);
  TH2F*  hCorrelationClustersTracklets = new TH2F("hCorrelationClustersTracklets","hCorrelationClustersTracklets;clusters ; tracklets ",500,0,500,100,0,100);
  TH2F*  hCorrelationHitsRadius = new TH2F("hCorrelationHitsRadius","hCorrelationHitsRadius;hits ; radius ",100,0,500,100,0,10);
  TH2F*  hCorrelationHitsX = new TH2F("hCorrelationHitsX","hCorrelationHitsX;hits ; X ",100,0,500,100,-5,5);
  TH2F*  hCorrelationHitsY = new TH2F("hCorrelationHitsY","hCorrelationHitsY;hits ; Y ",100,0,500,100,-5,5);
  fDiagList->Add(hCorrelationHitsRadius);
  fDiagList->Add(hCorrelationHitsX);
  fDiagList->Add(hCorrelationHitsY);
  fDiagList->Add(hCorrelationFMDSPD);
  fDiagList->Add(hCorrelationFMD1SPD);
  fDiagList->Add(hCorrelationFMD2ISPD);
  fDiagList->Add(hCorrelationFMD2OSPD);
  fDiagList->Add(hCorrelationFMD3ISPD);
  fDiagList->Add(hCorrelationFMD3OSPD);
  fDiagList->Add(hCorrelationFMDGoodtracks);
  fDiagList->Add(hCorrelationFMDBadtracks);
  fDiagList->Add(hCorrelationGoodbadtracks);
fDiagList->Add(hCorrelationFMDSPDhits);
  fDiagList->Add(hCorrelationClustersTracklets);
  fDiagList->Add(hCorrelationSPDTracklets);
  TH2F*  hCorrelationFMD23 = new TH2F("hCorrelationFMD23","hCorrelationFMD23;FMD2 ;FMD3 ",100,0,500,100,0,500);
  TH2F*  hCorrelationFMD2diff23 = new TH2F("hCorrelationFMD2diff23","hCorrelationFMD2diff23;FMD2 ;diff FMD23 ",100,0,100,100,0,100);
  TH2F*  hCorrelationFMD3diff23 = new TH2F("hCorrelationFMD3diff23","hCorrelationFMD3diff23;FMD3 ;diff FMD23 ",100,0,100,100,0,100);
  TH2F*  hCorrelationFMD1diff13 = new TH2F("hCorrelationFMD1diff13","hCorrelationFMD1diff13;FMD1 ;diff FMD13 ",100,0,100,100,0,100);
  TH2F*  hCorrelationFMD1diff12 = new TH2F("hCorrelationFMD1diff12","hCorrelationFMD1diff12;FMD1 ;diff FMD12 ",100,0,100,100,0,100);
  TH2F*  hCorrelationFMD12 = new TH2F("hCorrelationFMD12","hCorrelationFMD12;FMD1 ;FMD2 ",100,0,500,100,0,500);
  TH2F* hCorrelationFMDBgCand = new TH2F("hCorrelationFMDBgCand","hCorrelationFMDBgCand;Bg Tr ;FMD ",100,0,100,500,0,500);
  
  TH2F* hCorrelationFMDFlatTr = new TH2F("hCorrelationFMDFlatTr","hCorrelationFMDFlatTr;Bg Tr ;FMD ",100,0,100,500,0,500);
  TH2F* hCorrelationFMDRatioFlatTr = new TH2F("hCorrelationFMDRatioFlatTr","hCorrelationFMDRatioFlatTr;Bg Tr ;FMD ",100,0,1,500,0,500);
  fDiagList->Add(hCorrelationFMDBgCand);
  fDiagList->Add(hCorrelationFMDFlatTr);
  fDiagList->Add(hCorrelationFMDRatioFlatTr);
  TH2F* hCorrelationFMDBgCandRelative = new TH2F("hCorrelationFMDBgCandRelative","hCorrelationFMDBgCandRelative;Bg Tr ;FMD ",100,0,2,500,0,500);
  fDiagList->Add(hCorrelationFMDBgCandRelative);
  fDiagList->Add(hCorrelationFMD2diff23);
  fDiagList->Add(hCorrelationFMD3diff23);
  fDiagList->Add(hCorrelationFMD1diff13);
  fDiagList->Add(hCorrelationFMD1diff12);
  fDiagList->Add(hCorrelationFMD23);
  fDiagList->Add(hCorrelationFMD12);
  TH2F*  hTimeCorrelation = new TH2F("hCorrelationTime","hCorrelationTime ; time ; FMD hits",500,0,500,100,0,200);
  fDiagList->Add(hTimeCorrelation);
  TH1F*  hHitDistribution = new TH1F("hFMDHitDistribution","hFMDHitDistribution ; FMD hits",500,0,500);

  TH1F*  hHitDistributionFMD1 = new TH1F("hFMDHitDistributionFMD1","hFMDHitDistributionFMD1 ; FMD1 hits",500,0,500);
  TH1F*  hHitDistributionFMD2I = new TH1F("hFMDHitDistributionFMD2I","hFMDHitDistributionFMD2I ; FMD2I hits",500,0,500);
  TH1F*  hHitDistributionFMD2O = new TH1F("hFMDHitDistributionFMD2O","hFMDHitDistributionFMD2O ; FMD2O hits",500,0,500);
  TH1F*  hHitDistributionFMD3I = new TH1F("hFMDHitDistributionFMD3I","hFMDHitDistributionFMD3I ; FMD3I hits",500,0,500);
  TH1F*  hHitDistributionFMD3O = new TH1F("hFMDHitDistributionFMD3O","hFMDHitDistributionFMD3O ; FMD3O hits",500,0,500);
  TH1F*  hTrVtxDistribution = new TH1F("hTrVtxDistribution","hTrVtxDistribution ; TrVtx",200,-500,500);
  TH1F*  hTrEtaDistribution = new TH1F("hTrEtaDistribution","hTrEtaDistribution ; TrEta",200,-9,9);
  TH1F*  hTrEtaDistribution2 = new TH1F("hTrEtaDistribution2","hTrEtaDistribution2 ; TrEta",200,-9,9);
  
 TH1F*  hFlatTracks = new TH1F("hFlatTracks","hFlatTracks ; Horizontal tracks",100,0,100);
 
 TH1F* hEnergyOfParticles = new TH1F("hEnergyOfParticles","hEnergyOfParticles",1000000,0,10);
 fDiagList->Add(hEnergyOfParticles);
 fDiagList->Add(hTrVtxDistribution);
 fDiagList->Add(hTrEtaDistribution);
 fDiagList->Add(hTrEtaDistribution2);
  fDiagList->Add(hFlatTracks);
  fDiagList->Add(hHitDistribution);
  fDiagList->Add(hHitDistributionFMD1);
  fDiagList->Add(hHitDistributionFMD2I);
  fDiagList->Add(hHitDistributionFMD2O);
  fDiagList->Add(hHitDistributionFMD3I);
  fDiagList->Add(hHitDistributionFMD3O);
  TH1F*  nMCevents = new TH1F("nMCEventsNoCuts","nMCEventsNoCuts",pars->GetNvtxBins(),0,pars->GetNvtxBins());
  
  fDiagList->Add(nMCevents);
 
    
    
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::ConnectInputData(Option_t */*option*/)
{
  // connect the input data
  if(fStandalone)
    fESD = (AliESDEvent*)GetInputData(0);
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::Exec(Option_t */*option*/)
{
  //perform analysis on one event
  AliESD* old = fESD->GetAliESDOld();
  if (old) {
    fESD->CopyFromOldESD();
  }
  
  foutputESDFMD->Clear();
  
  Int_t delta = fESD->GetOrbitNumber() - fLastOrbit;
  fLastOrbit = fESD->GetOrbitNumber();
  
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  Double_t vertex[3];
  Bool_t vtxStatus = pars->GetVertex(fESD,vertex);
  fEsdVertex->SetXYZ(vertex);
  
  // Process primaries here to get true MC distribution
  if(pars->GetProcessPrimary())
    ProcessPrimary();
  const AliMultiplicity* testmult = fESD->GetMultiplicity();
  Int_t nTrackLets = testmult->GetNumberOfTracklets();
  TH2F*  hCorrelationClustersTracklets = (TH2F*)fDiagList->FindObject("hCorrelationClustersTracklets");
  hCorrelationClustersTracklets->Fill(testmult->GetNumberOfSingleClusters(),nTrackLets);  
  
  
  
  Bool_t isTriggered = pars->IsEventTriggered(fESD);
  
  if(!isTriggered || !vtxStatus ) {
    fStatus = kFALSE;
    return;
  }
  else
    fStatus = kTRUE;
  
  TH1F* hXvtx = (TH1F*)fDiagList->FindObject("hXvtx");
  if(vertex[0] != 0) hXvtx->Fill(vertex[0]);
  TH1F* hYvtx = (TH1F*)fDiagList->FindObject("hYvtx");
  if(vertex[1] != 0) hYvtx->Fill(vertex[1]);
  TH1F* hZvtx = (TH1F*)fDiagList->FindObject("hZvtx");
  hZvtx->Fill(vertex[2]);
  //const AliMultiplicity* testmult = fESD->GetMultiplicity();
  //std::cout<<vertex[2]<<std::endl;
  //Int_t nTrackLets = testmult->GetNumberOfTracklets();
  
  if( TMath::Abs(vertex[2]) > pars->GetVtxCutZ()) {
    fStatus = kFALSE;
    return;
  }
    
  if(nTrackLets < 1000) foutputESDFMD->SetUniqueID(kTRUE);
  else foutputESDFMD->SetUniqueID(kFALSE);
  
  AliESDFMD* fmd = fESD->GetFMDData();
  
  if (!fmd) return;
  Int_t nHits[3][2] = {{0,0},{0,0},{0,0}};

  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      
      TH1F* hEdist = (TH1F*)fDiagList->FindObject(Form("Edist_before_sharing_FMD%d%c",det,ring));
      
      for(UShort_t sec =0; sec < nsec;  sec++) {
	fSharedThis      = kFALSE;
	fSharedPrev      = kFALSE;
	
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  foutputESDFMD->SetMultiplicity(det,ring,sec,strip,0.);
	  Float_t mult = fmd->Multiplicity(det,ring,sec,strip);
	  
	  if(mult == AliESDFMD::kInvalidMult || mult == 0) continue;
	  
	 
	  //Double_t eta = fmd->Eta(det,ring,sec,strip);
	  Float_t eta = pars->GetEtaFromStrip(det,ring,sec,strip,vertex[2]);
	  	  
	  hEdist->Fill(mult);
	  if(fmd->IsAngleCorrected())
	    mult = mult/TMath::Cos(Eta2Theta(eta));
	  Float_t prevE = 0;
	  Float_t nextE = 0;
	  if(strip != 0)
	    if(fmd->Multiplicity(det,ring,sec,strip-1) != AliESDFMD::kInvalidMult) {
	      prevE = fmd->Multiplicity(det,ring,sec,strip-1);
	      if(fmd->IsAngleCorrected())
		prevE = prevE/TMath::Cos(Eta2Theta(fmd->Eta(det,ring,sec,strip-1)));
	    }
	  if(strip != nstr - 1)
	    if(fmd->Multiplicity(det,ring,sec,strip+1) != AliESDFMD::kInvalidMult) {
	      nextE = fmd->Multiplicity(det,ring,sec,strip+1);
	      if(fmd->IsAngleCorrected())
		nextE = nextE/TMath::Cos(Eta2Theta(fmd->Eta(det,ring,sec,strip+1)));
	    }
	  
	  Float_t mergedEnergy = GetMultiplicityOfStrip(mult,eta,prevE,nextE,det,ring,sec,strip);
	  //if(mult> (pars->GetMPV(det,ring,eta) - pars->GetSigma(det,ring,eta))) 
	  //  mergedEnergy = mult;
	  //else mergedEnergy = 0;
	  if(mergedEnergy > 0.3 )
	    nHits[det-1][ir]++;
	  foutputESDFMD->SetMultiplicity(det,ring,sec,strip,mergedEnergy);
	  foutputESDFMD->SetEta(det,ring,sec,strip,eta);
	  
	}
      }
    }
  }
  //cluster cut
  //if(testmult->GetNumberOfSingleClusters() > 15 + nTrackLets)
  //  {fStatus = kFALSE; std::cout<<"FMD : "<<nHits[0][0]<<"  "<<nHits[1][0]<<"   "<<nHits[1][1]<<"   "<<nHits[2][0]<<"   "<<nHits[2][1]<<" tracks  "<<testmult->GetNumberOfSingleClusters()<<"   "<<nTrackLets<<std::endl; return;}

  TH2F*  hCorrelationFMD23  = (TH2F*)fDiagList->FindObject("hCorrelationFMD23");
  TH2F*  hCorrelationFMD12  = (TH2F*)fDiagList->FindObject("hCorrelationFMD12");
  TH2F*  hCorrelationFMD2diff23  = (TH2F*)fDiagList->FindObject("hCorrelationFMD2diff23");
  TH2F*  hCorrelationFMD3diff23  = (TH2F*)fDiagList->FindObject("hCorrelationFMD3diff23");
  TH2F*  hCorrelationFMD1diff13  = (TH2F*)fDiagList->FindObject("hCorrelationFMD1diff13");  
  TH2F*  hCorrelationFMD1diff12  = (TH2F*)fDiagList->FindObject("hCorrelationFMD1diff12");

  TH2F*  hCorrelationFMDSPD = (TH2F*)fDiagList->FindObject("hCorrelationFMDSPD");
  TH2F*  hCorrelationFMD1SPD = (TH2F*)fDiagList->FindObject("hCorrelationFMD1SPD");
  TH2F*  hCorrelationFMD2ISPD = (TH2F*)fDiagList->FindObject("hCorrelationFMD2ISPD");
  TH2F*  hCorrelationFMD2OSPD = (TH2F*)fDiagList->FindObject("hCorrelationFMD2OSPD");
  TH2F*  hCorrelationFMD3ISPD = (TH2F*)fDiagList->FindObject("hCorrelationFMD3ISPD");
  TH2F*  hCorrelationFMD3OSPD = (TH2F*)fDiagList->FindObject("hCorrelationFMD3OSPD");
  
  TH2F*  hCorrelationFMDSPDhits = (TH2F*)fDiagList->FindObject("hCorrelationFMDSPDhits");

  TH2F*  hCorrelationSPDTracklets = (TH2F*)fDiagList->FindObject("hCorrelationSPDTracklets");
  TH2F*  hTimeCorrelation   = (TH2F*)fDiagList->FindObject("hCorrelationTime");
  TH2F*  hHitsRadius   = (TH2F*)fDiagList->FindObject("hCorrelationHitsRadius");
  TH2F*  hHitsX   = (TH2F*)fDiagList->FindObject("hCorrelationHitsX");
  TH2F*  hHitsY   = (TH2F*)fDiagList->FindObject("hCorrelationHitsY");
  TH1F*  hFMDHitDistribution   = (TH1F*)fDiagList->FindObject("hFMDHitDistribution");
  TH1F*  hFMD1HitDistribution   = (TH1F*)fDiagList->FindObject("hFMDHitDistributionFMD1");  
  TH1F*  hFMD2IHitDistribution   = (TH1F*)fDiagList->FindObject("hFMDHitDistributionFMD2I");
  TH1F*  hFMD2OHitDistribution   = (TH1F*)fDiagList->FindObject("hFMDHitDistributionFMD2O");
  TH1F*  hFMD3IHitDistribution   = (TH1F*)fDiagList->FindObject("hFMDHitDistributionFMD3I");
  TH1F*  hFMD3OHitDistribution   = (TH1F*)fDiagList->FindObject("hFMDHitDistributionFMD3O"); 
  TH1F*  hCorrelationFMDGoodtracks = (TH1F*)fDiagList->FindObject("hCorrelationFMDGoodtracks");
  TH1F*  hCorrelationFMDBadtracks = (TH1F*)fDiagList->FindObject("hCorrelationFMDBadtracks");
  TH1F*  hCorrelationGoodbadtracks = (TH1F*)fDiagList->FindObject("hCorrelationGoodbadtracks");
   TH2F*  hCorrelationFMDBgCand = (TH2F*)fDiagList->FindObject("hCorrelationFMDBgCand");
   TH2F*  hCorrelationFMDBgCandRelative = (TH2F*)fDiagList->FindObject("hCorrelationFMDBgCandRelative");
   
   TH2F* hCorrelationFMDFlatTr = (TH2F*)fDiagList->FindObject("hCorrelationFMDFlatTr");
   TH2F* hCorrelationFMDRatioFlatTr = (TH2F*)fDiagList->FindObject("hCorrelationFMDRatioFlatTr");
   TH1F* hTrVtxDistribution = (TH1F*)fDiagList->FindObject("hTrVtxDistribution");
   TH1F* hTrEtaDistribution = (TH1F*)fDiagList->FindObject("hTrEtaDistribution");
   TH1F* hTrEtaDistribution2 = (TH1F*)fDiagList->FindObject("hTrEtaDistribution2");
  hCorrelationFMDSPD->Fill(nTrackLets,nHits[0][0]+nHits[1][0]+nHits[1][1]+nHits[2][0]+nHits[2][1]);
  TH1F* hFlatTracks = (TH1F*)fDiagList->FindObject("hFlatTracks");
  hCorrelationFMD1SPD->Fill(nTrackLets,nHits[0][0]);
  hCorrelationFMD2ISPD->Fill(nTrackLets,nHits[1][0]);
  hCorrelationFMD2OSPD->Fill(nTrackLets,nHits[1][1]);
  hCorrelationFMD3ISPD->Fill(nTrackLets,nHits[2][0]);
  hCorrelationFMD3OSPD->Fill(nTrackLets,nHits[2][1]);
  hCorrelationFMDSPDhits->Fill(testmult->GetNumberOfFiredChips(0),nHits[0][0]+nHits[1][0]+nHits[1][1]+nHits[2][0]+nHits[2][1]);
  hCorrelationSPDTracklets->Fill(testmult->GetNumberOfFiredChips(0),nTrackLets);
  
  hTimeCorrelation->Fill(delta,nHits[0][0]+nHits[1][0]+nHits[1][1]+nHits[2][0]+nHits[2][1]);
  hCorrelationFMD23->Fill(nHits[1][0]+nHits[1][1],nHits[2][0]+nHits[2][1]);
  hCorrelationFMD12->Fill(nHits[0][0],nHits[1][0]+nHits[1][1]);
  
  //  if(TMath::Abs(nHits[1]-nHits[2]) > 15 && (nHits[1]+nHits[2]) > 35) {fStatus = kFALSE; std::cout<<"difference : "<<TMath::Abs(nHits[1]-nHits[2])<<std::endl; return;}
  
  // if(testmult->GetNumberOfFiredChips(0))
  //  if(testmult->GetNumberOfFiredChips(0) > 15 && ((Float_t)nTrackLets / (Float_t)testmult->GetNumberOfFiredChips(0)) < 0.4)
  //   {fStatus = kFALSE; std::cout<<nTrackLets<<"   "<<testmult->GetNumberOfFiredChips(0)<<"   "<<nHits[0]<<"   "<<nHits[1]<<"   "<<nHits[2]<<std::endl; return;}
  
  
  Float_t diff23 = (Float_t)TMath::Abs(nHits[2][1] + nHits[2][0] - nHits[1][1] - nHits[1][0]);
 
  Float_t diff13 = TMath::Abs(nHits[2][1] + nHits[2][0] - nHits[1][1] - nHits[1][0] - nHits[0][0]);
  Float_t diff12 = TMath::Abs(nHits[1][0] - nHits[0][0]);
  
  hCorrelationFMD1diff12->Fill(nHits[0][0], diff12);
  hCorrelationFMD1diff13->Fill(nHits[0][0], diff13);
  hCorrelationFMD2diff23->Fill(nHits[1][1], diff23);
  hCorrelationFMD3diff23->Fill(nHits[2][1], diff23);
  
  //
  Float_t nTotalFMDhits = nHits[0][0]+nHits[1][0]+nHits[1][1]+nHits[2][0]+nHits[2][1] ;
   Float_t radius = TMath::Sqrt(TMath::Power(vertex[0] + 0.03715,2) + TMath::Power(vertex[1] - 0.1659,2));
  
  if(vertex[1] !=0 || vertex[1] !=0) {
    hHitsRadius->Fill(nTotalFMDhits,radius);
    hHitsX->Fill(nTotalFMDhits,vertex[0]);
    hHitsY->Fill(nTotalFMDhits,vertex[1]); }
  
  hFMDHitDistribution->Fill(nTotalFMDhits);
  hFMD1HitDistribution->Fill(nHits[0][0]);
  hFMD2IHitDistribution->Fill(nHits[1][0]);
  hFMD2OHitDistribution->Fill(nHits[1][1]);
  hFMD3IHitDistribution->Fill(nHits[2][0]);
  hFMD3OHitDistribution->Fill(nHits[2][0]);
  
  // if(radius > 0.5) {fStatus = kFALSE; std::cout<<"FMD : "<<nTotalFMDhits<<std::endl; foutputESDFMD->Clear(); return;}
  //if(TMath::Abs(vertex[1] - 0.1659) > 0.1 ) {fStatus = kFALSE; std::cout<<"FMD : "<<nTotalFMDhits<<std::endl; foutputESDFMD->Clear(); return;}
  
   if(nTrackLets < pars->GetLowSPDLimit() || nTrackLets > pars->GetHighSPDLimit())
    {fStatus = kFALSE; std::cout<<nTrackLets<<"   "<<"   "<<nHits[0][0]<<"  "<<nHits[1][0]<<"   "<<nHits[1][1]<<"   "<<nHits[2][0]<<"   "<<nHits[2][1]<<std::endl; return;}
  
    
   AliESDtrack* track = 0;
  Int_t ntracks = fESD->GetNumberOfTracks();
  Float_t ngood =0, nbad = 0;
  //std::cout<<" Primary vtx : "<<vertex[0]<<"   "<<vertex[1]<<"   "<<vertex[2]<<"  "<<nTotalFMDhits<<std::endl;
  Int_t nBgCandidates = 0;
  Float_t nFlat = 0;
  for(Int_t i=0;i<ntracks;i++) {
    track = fESD->GetTrack(i);
    //std::cout<<track->GetX()-vertex[0]<<"   "<<track->GetY()-vertex[1]<<"   "<<track->GetZ()-vertex[2]<<std::endl;
    //std::cout<<track->GetX()<<"   "<<track->GetY()<<"   "<<track->GetZ()<<std::endl;
    hTrVtxDistribution->Fill(track->GetZ());
    
    if(TMath::Abs( track->GetZ()) > 50 &&  TMath::Abs(track->GetZ()) < 300) { // && TMath::Abs( track->GetY()) < 1)
      nBgCandidates++;
      hTrEtaDistribution->Fill(track->Eta()); 
    }
    
    
    
    if(TMath::Abs(track->GetX()-vertex[0]) > 0.3 || TMath::Abs(track->GetY()-vertex[1]) > 0.3  || TMath::Abs(track->GetZ()-vertex[2]) > 0.3)  {
      nbad++;
      hTrEtaDistribution2->Fill(track->Eta()); }
    else ngood++;
    
    if(TMath::Abs(track->Pt()) < 0.1)
      nFlat++;
  }
  
  Float_t ratioFlat = 0;
  if(fESD->GetNumberOfTracks())
    ratioFlat = nFlat/(Float_t)fESD->GetNumberOfTracks();
  hCorrelationFMDFlatTr->Fill(nFlat,nTotalFMDhits);
  hCorrelationFMDRatioFlatTr->Fill(ratioFlat,nTotalFMDhits);
  hFlatTracks->Fill(nFlat);
   
  // std::cout<<fESD->GetT0zVertex()<<"   "<<vertex[2]<<std::endl;
  Float_t ratioBg = 0;
  //if(fESD->GetNumberOfTracks() > 0)
  
  if(fESD->GetNumberOfTracks() > 0)
    ratioBg = (Float_t)nBgCandidates/(Float_t)fESD->GetNumberOfTracks();
  hCorrelationFMDBgCand->Fill(nBgCandidates,nTotalFMDhits);
  hCorrelationFMDBgCandRelative->Fill(ratioBg,nTotalFMDhits);
  
  
  // Float_t ratio =  (nbad > 0 ? ngood / nbad  : 0);
  
  hCorrelationFMDGoodtracks->Fill(ngood,nTotalFMDhits);
  hCorrelationFMDBadtracks->Fill(nbad,nTotalFMDhits);
  hCorrelationGoodbadtracks->Fill(ngood,nbad);
    
  if(fStandalone) {
    PostData(0, foutputESDFMD); 
    PostData(1, fEsdVertex); 
    PostData(2, fESD); 
    PostData(3, fDiagList); 
  }
}
//_____________________________________________________________________
Float_t AliFMDAnalysisTaskSharing::GetMultiplicityOfStrip(Float_t mult,
							  Float_t eta,
							  Float_t prevE,
							  Float_t nextE,
							  UShort_t   det,
							  Char_t  ring,
							  UShort_t /*sec*/,
							  UShort_t /*strip*/) {
  //analyse and perform sharing on one strip
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
 
  Float_t mergedEnergy = 0;
  //Float_t nParticles = 0;
  Float_t cutLow  = 0.3;//0.15;
  
  Float_t cutHigh = pars->GetMPV(det,ring,eta) - 2*pars->GetSigma(det,ring,eta);

  // if(ring == 'I')
  //  cutLow = 0.1;
  
  //cutLow = 0;
  //AliFMDParameters* recopars = AliFMDParameters::Instance();
  //cutLow = (5*recopars->GetPedestalWidth(det,ring,sec,strip))/(recopars->GetPulseGain(det,ring,sec,strip)*recopars->GetDACPerMIP());
  //if(foutputESDFMD->GetUniqueID() == kFALSE ) {
  
  if(mult > 12 || mult < cutLow)
    {
      //   std::cout<<"rejecting hit in FMD "<<det<<" "<<ring<<std::endl;
      fSharedThis      = kFALSE;
      fSharedPrev      = kFALSE;
      return 0;
    }
  
  
  
  
  // Float_t cutPart = pars->GetMPV(det,ring,eta) - 5*pars->GetSigma(det,ring,eta);
  Float_t totalE  = mult;
  

    //std::cout<<det<<ring<<"   "<<sec<<"    "<<strip<<"   "<<cutLow<<std::endl;
  if(fSharedThis) {
    fSharedThis      = kFALSE;
    fSharedPrev      = kTRUE;
    return 0.;
  }
  
  /*  if(mult < 0.33*pars->GetMPV(det,ring,eta)) {
    fSharedThis      = kFALSE;
    fSharedPrev      = kFALSE;
    return 0;
    }*/
  if(mult<nextE && nextE>cutHigh && foutputESDFMD->GetUniqueID() == kTRUE)
    {
      fSharedThis      = kFALSE;
      fSharedPrev      = kFALSE;
      return 0;
    }
  
  
  if(prevE > cutLow && prevE < cutHigh && !fSharedPrev ) {
    totalE += prevE;
  }
  
  if(nextE > cutLow && nextE < cutHigh ) {
    totalE += nextE;
    fSharedThis      = kTRUE;
  }
  totalE = totalE*TMath::Cos(Eta2Theta(eta));
  TH1F* hEdist = (TH1F*)fDiagList->FindObject(Form("Edist_after_sharing_FMD%d%c",det,ring));
  if(totalE > cutLow)
    hEdist->Fill(totalE);
  
  
  if(totalE > 0) {
    
    mergedEnergy = totalE;
    fSharedPrev      = kTRUE;
    // if(det == 1 && ring =='I')
  }
  else{
      fSharedThis      = kFALSE;
      fSharedPrev      = kFALSE;
  }

  // mergedEnergy = mult;
  

  /*   }
else {
    TH1F* hEdist = (TH1F*)fDiagList->FindObject(Form("Edist_after_sharing_FMD%d%c",det,ring));
    if(mult > cutLow)
      fEtotal+=mult;
    if(mult < cutLow) {
      mergedEnergy = fEtotal;
      fEtotal = 0;
      hEdist->Fill(mergedEnergy);
      
    }
    
   }*/
  
  return mergedEnergy; 
  //}  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::Terminate(Option_t* /* option*/) {
  
  TH1F*  hFMDHitDistribution   = (TH1F*)fDiagList->FindObject("hFMDHitDistribution");
  TH1F*  hFMD1HitDistribution   = (TH1F*)fDiagList->FindObject("hFMDHitDistributionFMD1");
  TH1F*  hFMD2IHitDistribution   = (TH1F*)fDiagList->FindObject("hFMDHitDistributionFMD2I");
  TH1F*  hFMD2OHitDistribution   = (TH1F*)fDiagList->FindObject("hFMDHitDistributionFMD2O");
  TH1F*  hFMD3IHitDistribution   = (TH1F*)fDiagList->FindObject("hFMDHitDistributionFMD3I");
  TH1F*  hFMD3OHitDistribution   = (TH1F*)fDiagList->FindObject("hFMDHitDistributionFMD3O");
  
  TH1F* hZvtx = (TH1F*)fDiagList->FindObject("hZvtx");
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      
      TH1F* hEdist = (TH1F*)fDiagList->FindObject(Form("Edist_before_sharing_FMD%d%c",det,ring));
      TH1F* hEdistAfter = (TH1F*)fDiagList->FindObject(Form("Edist_after_sharing_FMD%d%c",det,ring));
      if(hZvtx->GetEntries()) {
	hEdist->Scale(1./(Float_t)hZvtx->GetEntries());
	hEdistAfter->Scale(1./(Float_t)hZvtx->GetEntries());
      }
      
    }
    
  } 
  TH1F* hFlatTracks = (TH1F*)fDiagList->FindObject("hFlatTracks");
  TH1F* hTrVtxDistribution = (TH1F*)fDiagList->FindObject("hTrVtxDistribution");
  TH1F* hTrEtaDistribution = (TH1F*)fDiagList->FindObject("hTrEtaDistribution");
  TH1F* hTrEtaDistribution2 = (TH1F*)fDiagList->FindObject("hTrEtaDistribution2");
  if(hZvtx->GetEntries()) {
    hFMDHitDistribution->Scale(1./(Float_t)hZvtx->GetEntries());
    hFMD1HitDistribution->Scale(1./(Float_t)hZvtx->GetEntries());
    hFMD2IHitDistribution->Scale(1./(Float_t)hZvtx->GetEntries());
    hFMD2OHitDistribution->Scale(1./(Float_t)hZvtx->GetEntries());
    hFMD3IHitDistribution->Scale(1./(Float_t)hZvtx->GetEntries());
    hFMD3OHitDistribution->Scale(1./(Float_t)hZvtx->GetEntries());
    hFlatTracks->Scale(1./(Float_t)hZvtx->GetEntries());
    hTrVtxDistribution->Scale(1./(Float_t)hZvtx->GetEntries());
    hTrEtaDistribution->Scale(1./(Float_t)hZvtx->GetEntries());
    hTrEtaDistribution2->Scale(1./(Float_t)hZvtx->GetEntries());
  }
  

  
  
}
//_____________________________________________________________________
Float_t AliFMDAnalysisTaskSharing::Eta2Theta(Float_t eta) const{
  //convert the eta of a strip to a theta
  Float_t theta = 2*TMath::ATan(TMath::Exp(-1*eta));
  
  if(eta < 0)
    theta = theta-TMath::Pi();
  
  //  std::cout<<"From eta2Theta: "<<theta<<"   "<<eta<<std::endl;
  return theta;
  


}



//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::ProcessPrimary() {
  //Get the unspoiled MC dN/deta before event cuts
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  AliMCEvent* mcEvent = eventHandler->MCEvent();
  if(!mcEvent)
    return;
  fLastTrackByStrip.Reset(-1);
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  AliMCParticle* particle = 0;
  
  AliStack* stack = mcEvent->Stack();
  
  TH1F* hPrimary = (TH1F*)fDiagList->FindObject("hMultvsEtaNoCuts");
  TH1F* hEnergyOfParticles = (TH1F*)fDiagList->FindObject("hEnergyOfParticles");
  AliHeader* header            = mcEvent->Header();
  AliGenEventHeader* genHeader = header->GenEventHeader();
  
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
  
  if (!pythiaGenHeader) {
    std::cout<<" no pythia header!"<<std::endl;
    return; 
  }
  
	
  Int_t pythiaType = pythiaGenHeader->ProcessType();
  
  /*if(pythiaType==92||pythiaType==93){
      std::cout<<"single diffractive"<<std::endl;
      return;
     }
  if(pythiaType==94){
    std::cout<<"double diffractive"<<std::endl;
    return;
  }
  */
  std::cout<<pythiaType<<"   "<<stack->GetNprimary()<<std::endl;
  
  
  TArrayF vertex;
  genHeader->PrimaryVertex(vertex);
  
  if(TMath::Abs(vertex.At(2)) > pars->GetVtxCutZ())
    return;
  
  Double_t delta           = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
  Double_t vertexBinDouble = (vertex.At(2) + pars->GetVtxCutZ()) / delta;
  Int_t    vertexBin       = (Int_t)vertexBinDouble;
  
  Bool_t firstTrack = kTRUE;
  
  Int_t nTracks = stack->GetNprimary();
  if(pars->GetProcessHits())
    nTracks = stack->GetNtrack();
  TH1F* nMCevents = (TH1F*)fDiagList->FindObject("nMCEventsNoCuts");
  for(Int_t i = 0 ;i<nTracks;i++) {
    particle = (AliMCParticle*) mcEvent->GetTrack(i);
    if(!particle)
      continue;
    
    if(stack->IsPhysicalPrimary(i) && particle->Charge() != 0) {
      hPrimary->Fill(particle->Eta());
      

      TH1F* hPrimVtxBin = (TH1F*)fDiagList->FindObject(Form("primmult_NoCuts_vtxbin%d",vertexBin));
      hPrimVtxBin->Fill(particle->Eta());
      
      if(firstTrack) {
	nMCevents->Fill(vertexBin);
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
	if(particle->Charge() != 0)
	  hEnergyOfParticles->Fill(particle->E());
	
	AliFMDStripIndex::Unpack(ref->UserId(),det,ring,sec,strip);
	Float_t thisStripTrack = fLastTrackByStrip(det,ring,sec,strip);
	if(particle->Charge() != 0 && i != thisStripTrack ) {
	  //Double_t x,y,z;
	  
	  Float_t   eta   = pars->GetEtaFromStrip(det,ring,sec,strip,vertex.At(2));//-1*TMath::Log(TMath::Tan(0.5*theta));
	  TH1F* hHits = (TH1F*)fDiagList->FindObject(Form("hMCHits_nocuts_FMD%d%c_vtxbin%d",det,ring,vertexBin));
	  
	
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
// EOF
//
