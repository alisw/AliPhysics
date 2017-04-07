#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTriggerConfiguration.h"
#include "AliPID.h"
#include "AliCentrality.h"
#include "AliITSPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
#include <TParticle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include "AliESDInputHandlerRP.h"
#include "AliESDtrackCuts.h"
#include "AliMultSelection.h"
#include "AliAnalysisTaskITSsaTracks.h"

/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//*************************************************************************
// Implementation of class AliAnalysisTaskITSsaTracks
// AliAnalysisTaskSE to extract QA and performance histos for ITS standalone tracks
// 
//
// Authors: L. Milano, milano@to.infn.it
//          F. Prino, prino@to.infn.it
//          
//*************************************************************************

ClassImp(AliAnalysisTaskITSsaTracks)
//______________________________________________________________________________
AliAnalysisTaskITSsaTracks::AliAnalysisTaskITSsaTracks() : AliAnalysisTaskSE("ITSsa resolution"), 
  fOutput(0),
  fHistNEvents(0),
  fHistNEventsVsTrig(0),
  fHistPureSAtracksVsTracklets(0),
  fHistITSTPCtracksVsTracklets(0),
  fHistMCPhiResid(0),
  fHistPhiResid(0),
  fNtupleTracks(0),
  fNPtBins(kMaxPtBins),
  fMinITSpts(4),
  fMinSPDpts(1),
  fMinPtsforPid(3),
  fMinTPCpts(50),
  fMaxITSChi2Clu(100.),
  fMinCentrality(0.),
  fMaxCentrality(100.),
  fFillNtuple(kFALSE),
  fReadMC(kFALSE),
  fUseMCId(kFALSE),
  fUseCentrality(kFALSE),
  fRequireSPD(kTRUE),
  fRequireSDD(kTRUE),
  fRequireSSD(kTRUE),
  fTrigConfig(0),
  fInitCalib(kFALSE)
{
  //
  for(Int_t ilay=0; ilay<6;ilay++) fRequirePoint[ilay]=kFALSE;
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  const Int_t nbins = 35;
  Double_t xbins[nbins+1]={0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,
			   0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,
			   0.85,0.90,0.95,1.00,1.20,1.40,1.60,1.80,1.90,2.00,
			   3.00,4.00,5.00,10.0,20.0,30.0};
  SetPtBins(nbins,xbins);

  for(Int_t ij=0; ij<kNtrackTypes; ij++){
    fHistPt[ij]=0x0;
    fHistPtGood[ij]=0x0;
    fHistPtFake[ij]=0x0;
    for(Int_t ik=0; ik<7; ik++) fHistEtaPhiLay[ij*7+ik]=0x0;
    fHistEtaPhi[ij]=0x0;
    fHistEtaPhiGood[ij]=0x0;
    fHistEtaPhiFake[ij]=0x0;
    fHistEtaPhiAny[ij]=0x0;
    fHistEtaPhi1SPD[ij]=0x0;
    fHistEtaPhi4Clu[ij]=0x0;
    fHistEtaPhi6Clu[ij]=0x0;
    fHistChi2[ij]=0x0;
    fHistChi2Good[ij]=0x0;
    fHistChi2Fake[ij]=0x0;
    fHistNclu[ij]=0x0;
    fHistNcluGood[ij]=0x0;
    fHistNcluFake[ij]=0x0;
    fHistCluInLayVsPt[ij]=0x0;
    fHistdedxvsP2cls[ij]=0x0;
    fHistdedxvsP3cls[ij]=0x0;
    fHistdedxvsP4cls[ij]=0x0;
  }
 for(Int_t ij=0; ij<kNspecies; ij++){
   fHistPtTPCITS[ij]=0x0;
   fHistPtITSsa[ij]=0x0;
   fHistPtITSpureSA[ij]=0x0;
   fHistEtaPhiTPCITS[ij]=0x0;
   fHistEtaPhiITSsa[ij]=0x0;
   fHistEtaPhiITSpureSA[ij]=0x0;
   fHistNcluTPCITS[ij]=0x0;
   fHistNcluITSsa[ij]=0x0;
   fHistNcluITSpureSA[ij]=0x0;
   fHistd0rphiITSpureSA[ij]=0x0;
   fHistd0zITSpureSA[ij]=0x0;
   fHistCluInLayTPCITS[ij]=0x0;
   fHistCluInLayITSsa[ij]=0x0;
   fHistCluInLayITSpureSA[ij]=0x0;
   fHistOuterLayITSpureSA[ij]=0x0;
   fHistPtResid[ij]=0x0;
   fHistPtRelResid[ij]=0x0;
   fHistInvPtResid[ij]=0x0;
   fHistInvPtRelResid[ij]=0x0;
   fHistMCPtResid[ij]=0x0;
   fHistMCPtRelResid[ij]=0x0;
   fHistMCInvPtResid[ij]=0x0;
   fHistMCInvPtRelResid[ij]=0x0;
 }
}


//___________________________________________________________________________
AliAnalysisTaskITSsaTracks::~AliAnalysisTaskITSsaTracks(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistNEventsVsTrig;
    for(Int_t iType=0; iType<kNtrackTypes; iType++){
      delete fHistPt[iType];
      delete fHistPtGood[iType];
      delete fHistPtFake[iType];
      for(Int_t ik=0; ik<7; ik++) delete fHistEtaPhiLay[iType*7+ik];
      delete fHistEtaPhi[iType];
      delete fHistEtaPhiGood[iType];
      delete fHistEtaPhiFake[iType];
      delete fHistEtaPhiAny[iType];
      delete fHistEtaPhi1SPD[iType];
      delete fHistEtaPhi4Clu[iType];
      delete fHistEtaPhi6Clu[iType];
      delete fHistChi2[iType];
      delete fHistChi2Good[iType];
      delete fHistChi2Fake[iType];
      delete fHistNclu[iType];
      delete fHistNcluGood[iType];
      delete fHistNcluFake[iType];
      delete fHistCluInLayVsPt[iType];
      delete fHistdedxvsP2cls[iType];
      delete fHistdedxvsP3cls[iType];
      delete fHistdedxvsP4cls[iType];
    }
    for(Int_t iSpec=0; iSpec<kNspecies; iSpec++){
      delete fHistPtTPCITS[iSpec];
      delete fHistPtITSsa[iSpec];
      delete fHistPtITSpureSA[iSpec];
      delete fHistEtaPhiTPCITS[iSpec];
      delete fHistEtaPhiITSsa[iSpec];
      delete fHistEtaPhiITSpureSA[iSpec];
      delete fHistNcluTPCITS[iSpec];
      delete fHistNcluITSsa[iSpec];
      delete fHistNcluITSpureSA[iSpec];
      delete fHistd0rphiITSpureSA[iSpec];
      delete fHistd0zITSpureSA[iSpec];
      delete fHistCluInLayTPCITS[iSpec];
      delete fHistCluInLayITSsa[iSpec];
      delete fHistCluInLayITSpureSA[iSpec];
      delete fHistOuterLayITSpureSA[iSpec];
      delete fHistPtResid[iSpec];
      delete fHistPtRelResid[iSpec];
      delete fHistInvPtResid[iSpec];
      delete fHistInvPtRelResid[iSpec];
      delete fHistMCPtResid[iSpec];
      delete fHistMCPtRelResid[iSpec];
      delete fHistMCInvPtResid[iSpec];
      delete fHistMCInvPtRelResid[iSpec];
    } 
  delete fHistMCPhiResid;
  delete fHistPhiResid;
  delete fNtupleTracks;
  }
  delete fOutput;
}
   
//________________________________________________________________________
void AliAnalysisTaskITSsaTracks::SetPtBins(Int_t n, Double_t* lim){
  // define pt bins for analysis
  if(n>kMaxPtBins){
    printf("Max. number of Pt bins = %d\n",kMaxPtBins);
    return;
  }else{
    fNPtBins=n;
    for(Int_t i=0; i<fNPtBins+1; i++) fPtLimits[i]=lim[i];
    for(Int_t i=fNPtBins+1; i<kMaxPtBins+1; i++) fPtLimits[i]=99999999.;
  }
}
//___________________________________________________________________________
void AliAnalysisTaskITSsaTracks::UserCreateOutputObjects() {
  // create output histos

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fNtupleTracks= new TNtuple("ntupleTracks","ntupleTracks","pt:px:py:pz:eta:phi:d0xy:d0z:dedx3:dedx4:dedx5:dedx6:truncmean:chi2:status:ITSrefit:TPCin:TPCrefit:ITSpureSA:nITSclu:nTPCclu:clumap:spdin:spdout:sddin:sddout:ssdin:ssdout:label:ptgen:pxgen:pygen:pzgen:etagen:phigen:ntracks:ntracklets:centr:xvert:yvert:zvert:iEvent");
  // kinematics: pt, p eta,phi        ->  6 variables
  // impact parameters: d0xy, d0z     ->  2 variables
  // dE/dx: 4 Layers + trunc mean     ->  5 variables
  // chi2 and track status:           ->  6 variables
  // cluster infos:                   ->  9 variables
  // MC info:                         ->  7 variables
  // Multiplicity                     ->  1 variable
  // Total:                              36
  fOutput->Add(fNtupleTracks);

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",15,-0.5,14.5);
  //fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"All events");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"Good vertex"); 
  fHistNEvents->GetXaxis()->SetBinLabel(3,"In centrality range"); 
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Without SPD"); 
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Without SDD"); 
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Without SSD"); 
  fHistNEvents->GetXaxis()->SetBinLabel(7,"No ITS layers"); 
  fHistNEvents->GetXaxis()->SetBinLabel(8,"SPD only"); 
  fHistNEvents->GetXaxis()->SetBinLabel(9,"SDD only"); 
  fHistNEvents->GetXaxis()->SetBinLabel(10,"SSD only"); 
  fHistNEvents->GetXaxis()->SetBinLabel(11,"SPD+SDD");
  fHistNEvents->GetXaxis()->SetBinLabel(12,"SPD+SSD");
  fHistNEvents->GetXaxis()->SetBinLabel(13,"SDD+SSD");
  fHistNEvents->GetXaxis()->SetBinLabel(14,"SPD+SDD+SSD");
  fHistNEvents->GetXaxis()->SetBinLabel(15,"Analyzed events"); 
  fOutput->Add(fHistNEvents);

  fHistNEventsVsTrig = new TH2F("hNEventsVsTrig", "Events per trigger type",9,-0.5,8.5,8,-0.5,7.5);
  fHistNEventsVsTrig->GetXaxis()->SetBinLabel(1,"MB - all");
  fHistNEventsVsTrig->GetXaxis()->SetBinLabel(2,"MB - fast");
  fHistNEventsVsTrig->GetXaxis()->SetBinLabel(3,"HighMult - all");
  fHistNEventsVsTrig->GetXaxis()->SetBinLabel(4,"HighMult - fast");
  fHistNEventsVsTrig->GetXaxis()->SetBinLabel(5,"Central");
  fHistNEventsVsTrig->GetXaxis()->SetBinLabel(6,"SemiCentral");
  fHistNEventsVsTrig->GetXaxis()->SetBinLabel(7,"EMC");
  fHistNEventsVsTrig->GetXaxis()->SetBinLabel(8,"Muon");
  fHistNEventsVsTrig->GetXaxis()->SetBinLabel(9,"Other");
  fHistNEventsVsTrig->GetYaxis()->SetBinLabel(1,"No ITS layers"); 
  fHistNEventsVsTrig->GetYaxis()->SetBinLabel(2,"SPD only"); 
  fHistNEventsVsTrig->GetYaxis()->SetBinLabel(3,"SDD only"); 
  fHistNEventsVsTrig->GetYaxis()->SetBinLabel(4,"SSD only"); 
  fHistNEventsVsTrig->GetYaxis()->SetBinLabel(5,"SPD+SDD");
  fHistNEventsVsTrig->GetYaxis()->SetBinLabel(6,"SPD+SSD");
  fHistNEventsVsTrig->GetYaxis()->SetBinLabel(7,"SDD+SSD");
  fHistNEventsVsTrig->GetYaxis()->SetBinLabel(8,"SPD+SDD+SSD");
  fOutput->Add(fHistNEventsVsTrig);


  fHistPureSAtracksVsTracklets = new TH2F("hPureSAtracksVsTracklets","  ; nTrackelts ; nITSpureSAtracks",100,0.,200,100,0.,200.);
  fOutput->Add(fHistPureSAtracksVsTracklets);

  fHistITSTPCtracksVsTracklets = new TH2F("hITSTPCtracksVsTracklets","  ; nTrackelts ; nITSTPCtracks",100,0.,200,100,0.,200.);
  fOutput->Add(fHistITSTPCtracksVsTracklets);

  //binning for the dedx histogram
  const Int_t hnbinsdedx=400;
  Double_t hxmindedx = 0.01;
  Double_t hxmaxdedx = 10;
  Double_t hlogxmindedx = TMath::Log10(hxmindedx);
  Double_t hlogxmaxdedx = TMath::Log10(hxmaxdedx);
  Double_t hbinwidthdedx = (hlogxmaxdedx-hlogxmindedx)/hnbinsdedx;
  Double_t hxbinsdedx[hnbinsdedx+1];
  hxbinsdedx[0] = 0.01; 
  for (Int_t i=1;i<=hnbinsdedx;i++) {
    hxbinsdedx[i] = hxmindedx + TMath::Power(10,hlogxmindedx+i*hbinwidthdedx);
  }
  

  TString tit[kNtrackTypes]={"TPCITS","ITSsa","ITSpureSA"};
  for(Int_t iType=0; iType<kNtrackTypes; iType++){
    fHistPt[iType]=new TH1F(Form("hPt%s",tit[iType].Data()),"",100,0.,2.);
    //fHistPt[iType]->Sumw2();
    fOutput->Add(fHistPt[iType]);
    fHistPtGood[iType]=new TH1F(Form("hPtGood%s",tit[iType].Data()),"",100,0.,2.);
    //fHistPtGood[iType]->Sumw2();
    fOutput->Add(fHistPtGood[iType]);
    fHistPtFake[iType]=new TH1F(Form("hPtFake%s",tit[iType].Data()),"",100,0.,2.);
    //fHistPtFake[iType]->Sumw2();
    fOutput->Add(fHistPtFake[iType]);

    fHistEtaPhiLay[7*iType]=new TH2F(Form("hEtaPhiTracksNoLaySel%s",tit[iType].Data()),"",50,-1.,1.,200,0.,2.*TMath::Pi());
    fHistEtaPhiLay[7*iType]->SetMinimum(0);
    fOutput->Add(fHistEtaPhiLay[7*iType]);
    for(Int_t iLay=1; iLay<=6; iLay++){
      fHistEtaPhiLay[7*iType+iLay]=new TH2F(Form("hEtaPhiTracksLay%d%s",iLay,tit[iType].Data()),"",50,-1.,1.,200,0.,2.*TMath::Pi());
      fHistEtaPhiLay[7*iType+iLay]->SetMinimum(0);
      fOutput->Add(fHistEtaPhiLay[7*iType+iLay]);
    }
    
    fHistEtaPhi[iType] = new TH2F(Form("hEtaPhi%s",tit[iType].Data()),"",50,-1,1.,100,0.,2.*TMath::Pi());
    //fHistEtaPhi[iType]->Sumw2();
    fOutput->Add(fHistEtaPhi[iType]);
    fHistEtaPhiGood[iType] = new TH2F(Form("hEtaPhiGood%s",tit[iType].Data()),"",50,-1,1.,100,0.,2.*TMath::Pi());
    //fHistEtaPhiGood[iType]->Sumw2();
    fOutput->Add(fHistEtaPhiGood[iType]);
    fHistEtaPhiFake[iType] = new TH2F(Form("hEtaPhiFake%s",tit[iType].Data()),"",50,-1,1.,100,0.,2.*TMath::Pi());
    //fHistEtaPhiFake[iType]->Sumw2();
    fOutput->Add(fHistEtaPhiFake[iType]);

    fHistEtaPhiAny[iType] = new TH2F(Form("hEtaPhiAny%s",tit[iType].Data()),"",50,-1,1.,100,0.,2.*TMath::Pi());
    //fHistEtaPhiAny[iType]->Sumw2();
    fOutput->Add(fHistEtaPhiAny[iType]);
    fHistEtaPhi1SPD[iType] = new TH2F(Form("hEtaPhi1SPD%s",tit[iType].Data()),"",50,-1,1.,100,0.,2.*TMath::Pi());
    //fHistEtaPhi1SPD[iType]->Sumw2();
    fOutput->Add(fHistEtaPhi1SPD[iType]);
    fHistEtaPhi4Clu[iType] = new TH2F(Form("hEtaPhi4Clu%s",tit[iType].Data()),"",50,-1,1.,100,0.,2.*TMath::Pi());
    //fHistEtaPhi4Clu[iType]->Sumw2();
    fOutput->Add(fHistEtaPhi4Clu[iType]);
    fHistEtaPhi6Clu[iType] = new TH2F(Form("hEtaPhi6Clu%s",tit[iType].Data()),"",50,-1,1.,100,0.,2.*TMath::Pi());
    //fHistEtaPhi6Clu[iType]->Sumw2();
    fOutput->Add(fHistEtaPhi6Clu[iType]);

    fHistChi2[iType]=new TH1F(Form("hChi2%s",tit[iType].Data()),"",100,0.,10.);
    //fHistChi2[iType]->Sumw2();
    fOutput->Add(fHistChi2[iType]);
    fHistChi2Good[iType]=new TH1F(Form("hChi2Good%s",tit[iType].Data()),"",100,0.,10.);
    //fHistChi2Good[iType]->Sumw2();
    fOutput->Add(fHistChi2Good[iType]);
    fHistChi2Fake[iType]=new TH1F(Form("hChi2Fake%s",tit[iType].Data()),"",100,0.,10.);
    //fHistChi2Fake[iType]->Sumw2();
    fOutput->Add(fHistChi2Fake[iType]);

    fHistNclu[iType]=new TH1F(Form("hNclu%s",tit[iType].Data()),"",7,-0.5,6.5);
    //fHistNclu[iType]->Sumw2();
    fOutput->Add(fHistNclu[iType]);
    fHistNcluGood[iType]=new TH1F(Form("hNcluGood%s",tit[iType].Data()),"",7,-0.5,6.5);
    //fHistNcluGood[iType]->Sumw2();
    fOutput->Add(fHistNcluGood[iType]);
    fHistNcluFake[iType]=new TH1F(Form("hNcluFake%s",tit[iType].Data()),"",7,-0.5,6.5);
    //fHistNcluFake[iType]->Sumw2();
    fOutput->Add(fHistNcluFake[iType]);

    fHistCluInLayVsPt[iType]=new TH2F(Form("hCluInLayVsPt%s",tit[iType].Data()),"",100,0.,2.,7,-1.5,5.5);
    fOutput->Add(fHistCluInLayVsPt[iType]);

    fHistdedxvsP2cls[iType] = new TH2F(Form("hdedxvsP2cls%s",tit[iType].Data()),"",hnbinsdedx,hxbinsdedx,900,0,1000);
    //fHistdedxvsP2cls[iType]->Sumw2();
    fOutput->Add(fHistdedxvsP2cls[iType]);

    fHistdedxvsP3cls[iType] = new TH2F(Form("hdedxvsP3cls%s",tit[iType].Data()),"",hnbinsdedx,hxbinsdedx,900,0,1000);
    //fHistdedxvsP3cls[iType]->Sumw2();
    fOutput->Add(fHistdedxvsP3cls[iType]);

    fHistdedxvsP4cls[iType] = new TH2F(Form("hdedxvsP4cls%s",tit[iType].Data()),"",hnbinsdedx,hxbinsdedx,900,0,1000);
    //fHistdedxvsP4cls[iType]->Sumw2();
    fOutput->Add(fHistdedxvsP4cls[iType]);
  }

  //-----------------------------------------------------------

  TString spname[3]={"Pion","Kaon","Proton"};
  TString hisname;
  const Int_t nbins = fNPtBins;
  Double_t xbins[nbins+1];
  for(Int_t ibin=0; ibin<=nbins; ibin++) xbins[ibin]=fPtLimits[ibin];

  for(Int_t iSpec=0; iSpec<kNspecies; iSpec++){

    hisname.Form("hPtTPCITS%s",spname[iSpec].Data());
    fHistPtTPCITS[iSpec] = new TH1F(hisname.Data(),"",100,0.,2.);
    //fHistPtTPCITS[iSpec]->Sumw2();
    fOutput->Add(fHistPtTPCITS[iSpec]);

    hisname.Form("hPtITSsa%s",spname[iSpec].Data());
    fHistPtITSsa[iSpec] = new TH1F(hisname.Data(),"",100,0.,2.);
    //fHistPtITSsa[iSpec]->Sumw2();
    fOutput->Add(fHistPtITSsa[iSpec]);

    hisname.Form("hPtITSpureSA%s",spname[iSpec].Data());
    fHistPtITSpureSA[iSpec] = new TH1F(hisname.Data(),"",100,0.,2.);
    //fHistPtITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistPtITSpureSA[iSpec]);

    //---

    hisname.Form("hEtaPhiTPCITS%s",spname[iSpec].Data());
    fHistEtaPhiTPCITS[iSpec] = new TH2F(hisname.Data(),"",50,-1,1.,50,0.,2.*TMath::Pi());
    //fHistEtaPhiTPCITS[iSpec]->Sumw2();
    fOutput->Add(fHistEtaPhiTPCITS[iSpec]);

    hisname.Form("hEtaPhiITSsa%s",spname[iSpec].Data());
    fHistEtaPhiITSsa[iSpec] = new TH2F(hisname.Data(),"",50,-1,1.,50,0.,2.*TMath::Pi());
    //fHistEtaPhiITSsa[iSpec]->Sumw2();
    fOutput->Add(fHistEtaPhiITSsa[iSpec]);

    hisname.Form("hEtaPhiITSpureSA%s",spname[iSpec].Data());
    fHistEtaPhiITSpureSA[iSpec] = new TH2F(hisname.Data(),"",50,-1,1.,50,0.,2.*TMath::Pi());
    //fHistEtaPhiITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistEtaPhiITSpureSA[iSpec]);

    //---

    hisname.Form("hNcluTPCITS%s",spname[iSpec].Data());
    fHistNcluTPCITS[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-0.5,6.5);
    //fHistNcluTPCITS[iSpec]->Sumw2();
    fOutput->Add(fHistNcluTPCITS[iSpec]);

    hisname.Form("hNcluITSsa%s",spname[iSpec].Data());
    fHistNcluITSsa[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-0.5,6.5);
    //fHistNcluITSsa[iSpec]->Sumw2();
    fOutput->Add(fHistNcluITSsa[iSpec]);

    hisname.Form("hNcluITSpureSA%s",spname[iSpec].Data());
    fHistNcluITSpureSA[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-0.5,6.5);
    //fHistNcluITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistNcluITSpureSA[iSpec]);

    //---

    hisname.Form("hd0rphiITSpureSA%s",spname[iSpec].Data());
    fHistd0rphiITSpureSA[iSpec] = new TH2F(hisname.Data(),"",nbins,xbins,4000,-2,2);
    //fHistd0rphiITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistd0rphiITSpureSA[iSpec]);

    hisname.Form("hd0zITSpureSA%s",spname[iSpec].Data());
    fHistd0zITSpureSA[iSpec] = new TH2F(hisname.Data(),"",nbins,xbins,4000,-2,2);
    //fHistd0zITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistd0zITSpureSA[iSpec]);

    //---

    hisname.Form("hCluInLayTPCITS%s",spname[iSpec].Data());
    fHistCluInLayTPCITS[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-1.5,5.5);
    //fHistCluInLayTPCITS[iSpec]->Sumw2();
    fOutput->Add(fHistCluInLayTPCITS[iSpec]);
    
    hisname.Form("hCluInLayITSsa%s",spname[iSpec].Data());
    fHistCluInLayITSsa[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-1.5,5.5);
    //fHistCluInLayITSsa[iSpec]->Sumw2();
    fOutput->Add(fHistCluInLayITSsa[iSpec]);

    hisname.Form("hCluInLayITSpureSA%s",spname[iSpec].Data());
    fHistCluInLayITSpureSA[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-1.5,5.5);
    //fHistCluInLayITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistCluInLayITSpureSA[iSpec]);
    
    hisname.Form("hOuterLayITSpureSA%s",spname[iSpec].Data());
    fHistOuterLayITSpureSA[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-1.5,5.5);
    //fHistOuterLayITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistOuterLayITSpureSA[iSpec]);
    
    //---

    hisname.Form("hPtResid%s",spname[iSpec].Data());
    fHistPtResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-1.,1.);
    //fHistPtResid[iSpec]->Sumw2();
    fOutput->Add(fHistPtResid[iSpec]);
    hisname.Form("hPtRelResid%s",spname[iSpec].Data());
    fHistPtRelResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-0.5,0.5);
    //fHistPtRelResid[iSpec]->Sumw2();
    fOutput->Add(fHistPtRelResid[iSpec]);
    
    hisname.Form("hInvPtResid%s",spname[iSpec].Data());
    fHistInvPtResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-1.,1.);
    //fHistInvPtResid[iSpec]->Sumw2();
    fOutput->Add(fHistInvPtResid[iSpec]);
    hisname.Form("hInvPtRelResid%s",spname[iSpec].Data());
    fHistInvPtRelResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-0.5,0.5);
    //fHistInvPtRelResid[iSpec]->Sumw2();
    fOutput->Add(fHistInvPtRelResid[iSpec]);
    
    hisname.Form("hMCPtResid%s",spname[iSpec].Data());
    fHistMCPtResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-1.,1.);
    //fHistMCPtResid[iSpec]->Sumw2();
    fOutput->Add(fHistMCPtResid[iSpec]);
    hisname.Form("hMCPtRelResid%s",spname[iSpec].Data());
    fHistMCPtRelResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-0.5,0.5);
    //fHistMCPtRelResid[iSpec]->Sumw2();
    fOutput->Add(fHistMCPtRelResid[iSpec]);
    
    hisname.Form("hMCInvPtResid%s",spname[iSpec].Data());
    fHistMCInvPtResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-1.,1.);
    //fHistMCInvPtResid[iSpec]->Sumw2();
    fOutput->Add(fHistMCInvPtResid[iSpec]);
    hisname.Form("hMCInvPtRelResid%s",spname[iSpec].Data());
    fHistMCInvPtRelResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-0.5,0.5);
    //fHistMCInvPtRelResid[iSpec]->Sumw2();
    fOutput->Add(fHistMCInvPtRelResid[iSpec]);
    
  }

  fHistMCPhiResid=new TH2F("hMCPhiResid","",nbins,xbins,100,-0.5,0.5);
  //fHistMCPhiResid->Sumw2();
  fOutput->Add(fHistMCPhiResid);
  fHistPhiResid=new TH2F("hPhiResid","",nbins,xbins,100,-0.5,0.5);
  //fHistPhiResid->Sumw2();
  fOutput->Add(fHistPhiResid);

  PostData(1,fOutput);

}
//______________________________________________________________________________
void AliAnalysisTaskITSsaTracks::UserExec(Option_t *)
{
  //

  AliESDEvent *esd = (AliESDEvent*) (InputEvent());
  if(!esd) {
    printf("AliAnalysisTaskITSsaTracks::UserExec(): bad ESD\n");
    return;
  } 

  if(!ESDfriend()) {
    printf("AliAnalysisTaskITSsaTracks::UserExec(): bad ESDfriend\n");
    return;
  }

  if(fRequireSPD || fRequireSDD || fRequireSSD){
    if (!fInitCalib) {
      AliCDBManager* man = AliCDBManager::Instance();
      if (!man) {
	AliFatal("CDB not set but needed by AliAnalysisTaskITSsaTracks");
	return;
      }   
      AliCDBEntry* eT=(AliCDBEntry*)man->Get("GRP/CTP/Config");
      if(eT){
	eT->PrintId();
	eT->PrintMetaData();
	fTrigConfig=(AliTriggerConfiguration*)eT->GetObject();
      }
      if(!eT || !fTrigConfig){
	AliError("Cannot retrieve CDB entry for GRP/CTP/Config");
	return;      
      }
      fInitCalib=kTRUE;
    }
  }

  AliMCEvent* mcEvent = 0x0;

  if(fReadMC){
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }
    mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
  }


  fHistNEvents->Fill(0);
  const AliESDVertex *spdv=esd->GetPrimaryVertexSPD();
  if(spdv->GetNContributors()<=0) return;
  fHistNEvents->Fill(1);

  Float_t xvert=spdv->GetX();
  Float_t yvert=spdv->GetY();
  Float_t zvert=spdv->GetZ();

  const Int_t ntSize=42;
  Float_t xnt[ntSize];
  
  Double_t centrality=-1.;
  if(esd->GetRunNumber()<244824){
    centrality=esd->GetCentrality()->GetCentralityPercentile("V0M");
  }else{
    AliMultSelection *multSelection = (AliMultSelection * )esd->FindListObject("MultSelection");
    if(multSelection){
      centrality=multSelection->GetMultiplicityPercentile("V0M"); 
      Int_t qual = multSelection->GetEvSelCode();
      if(qual == 199 ) centrality=-1.;
    }
  }

  if(fUseCentrality){
    if(centrality<fMinCentrality) return;
    if(centrality>fMaxCentrality) return;
  }
  fHistNEvents->Fill(2);

  Int_t trigMask=esd->GetTriggerMask();
  Bool_t isMB=trigMask&(AliVEvent::kAnyINT);
  Int_t muonMask=(AliVEvent::kMUON) | (AliVEvent::kCMUS5) | (AliVEvent::kMUSH7) | (AliVEvent::kMUL7) | (AliVEvent::kMUU7) | (AliVEvent::kMUS7) | (AliVEvent::kMuonSingleLowPt8) | (AliVEvent::kMuonSingleHighPt8) | (AliVEvent::kMuonLikeLowPt8) | (AliVEvent::kMuonUnlikeLowPt8) | (AliVEvent::kMuonUnlikeLowPt0);
  Bool_t isMuon=trigMask&muonMask;
  Bool_t isHighMult=trigMask&(AliVEvent::kHighMult);
  Int_t emcMask=(AliVEvent::kEMC1) | (AliVEvent::kEMC7) | (AliVEvent::kEMC8) | (AliVEvent::kEMCEJE) | (AliVEvent::kEMCEGA);
  Bool_t isEMC=trigMask&emcMask;
  Bool_t isCent=trigMask&(AliVEvent::kCentral);
  Bool_t isSemiCent=trigMask&(AliVEvent::kSemiCentral);
  Bool_t isFast=trigMask&(AliVEvent::kFastOnly);
  Bool_t spdOK=kTRUE;
  Bool_t sddOK=kTRUE;
  Bool_t ssdOK=kTRUE;
  if(fTrigConfig){
    spdOK=esd->IsDetectorInTriggerCluster("ITSSPD",fTrigConfig);
    if(!spdOK) fHistNEvents->Fill(3);
    sddOK=esd->IsDetectorInTriggerCluster("ITSSDD",fTrigConfig);
    if(!sddOK) fHistNEvents->Fill(4);
    ssdOK=esd->IsDetectorInTriggerCluster("ITSSSD",fTrigConfig);
    if(!ssdOK) fHistNEvents->Fill(5);
    Int_t laymask=0;
    if(spdOK) laymask+=1;
    if(sddOK) laymask+=2;
    if(ssdOK) laymask+=4;
    fHistNEvents->Fill(6+laymask);
    if(isMB){
      if(!isFast) fHistNEventsVsTrig->Fill(0.,laymask);
      else fHistNEventsVsTrig->Fill(1.,laymask);
    }
    else if(isHighMult){
      if(!isFast) fHistNEventsVsTrig->Fill(2.,laymask);
      else fHistNEventsVsTrig->Fill(3.,laymask);
    }
    else if(isCent) fHistNEventsVsTrig->Fill(4.,laymask);
    else if(isSemiCent) fHistNEventsVsTrig->Fill(5.,laymask);
    else if(isEMC) fHistNEventsVsTrig->Fill(6.,laymask);
    else if(isMuon) fHistNEventsVsTrig->Fill(7.,laymask);
    else fHistNEventsVsTrig->Fill(8.,laymask);
  }
  if(fRequireSPD && !spdOK) return;
  if(fRequireSDD && !sddOK) return;
  if(fRequireSSD && !ssdOK) return;
  fHistNEvents->Fill(14);
  

  Int_t ntracks = esd->GetNumberOfTracks();
  const AliMultiplicity* mult=esd->GetMultiplicity();
  Int_t ntracklets = mult->GetNumberOfTracklets();
  Int_t nPureSAtracks=0;
  Int_t nITSTPCtracks=0;

  AliESDtrackCuts* esdTrackCutsTPC = new AliESDtrackCuts("esdtrackCutsTPC");
  esdTrackCutsTPC->SetMinNCrossedRowsTPC(50);
  esdTrackCutsTPC->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  // esdTrackCutsTPC->SetCutGeoNcrNcl(2., 130., 1.5, 0.0, 0.0); 
  // esdTrackCutsTPC->SetCutOutDistortedRegionsTPC(kTRUE);
  esdTrackCutsTPC->SetMaxChi2PerClusterTPC(4);
  esdTrackCutsTPC->SetAcceptKinkDaughters(kFALSE);
  esdTrackCutsTPC->SetRequireTPCRefit(kTRUE);
  esdTrackCutsTPC->SetDCAToVertex2D(kFALSE);
  esdTrackCutsTPC->SetRequireSigmaToVertex(kFALSE);
  esdTrackCutsTPC->SetMaxChi2TPCConstrainedGlobal(36);
  esdTrackCutsTPC->SetMaxChi2PerClusterITS(36);

  for (Int_t iTrack=0; iTrack < ntracks; iTrack++) {
    AliESDtrack * track = esd->GetTrack(iTrack);
    if (!track) continue;
    Float_t pttrack=track->Pt();
    Float_t ptrack=track->P();
    Float_t pxtrack=track->Px();
    Float_t pytrack=track->Py();
    Float_t pztrack=track->Pz();
    Float_t impactXY=-999, impactZ=-999;
    track->GetImpactParameters(impactXY, impactZ);
    Double_t dedxlay[4];
    track->GetITSdEdxSamples(dedxlay);
    Float_t dedx=track->GetITSsignal();
    Float_t chi2=track->GetITSchi2();
    Int_t status=track->GetStatus();
    Bool_t isITSrefit=kFALSE;
    Bool_t isTPCin=kFALSE;
    Bool_t isTPCrefit=kFALSE;
    Bool_t isPureSA=kFALSE;
    if(status&AliESDtrack::kITSrefit) isITSrefit=kTRUE; 
    if(status&AliESDtrack::kTPCin){ 
      isTPCin=kTRUE; 
      if(!esdTrackCutsTPC->AcceptTrack(track)) continue;
    }
    if(status&AliESDtrack::kTPCrefit) isTPCrefit=kTRUE; 
    if(status&AliESDtrack::kITSpureSA) isPureSA=kTRUE;
    Bool_t isSA=kTRUE;
    if(status&AliESDtrack::kTPCin) isSA=kFALSE;
    Int_t nTPCclus=track->GetNcls(1);
    Int_t nITSclus=track->GetNcls(0);
    UChar_t clumap=track->GetITSClusterMap();
    Int_t idet;
    Float_t xloc,zloc;
    Int_t statusLay[6];
    for(Int_t iLay=0; iLay<6;iLay++) track->GetITSModuleIndexInfo(iLay,idet,statusLay[iLay],xloc,zloc);
    Int_t trlabel=track->GetLabel();
    Float_t ptgen=-999.;
    //    Float_t pgen=-999.;
    Float_t pxgen=-999.;
    Float_t pygen=-999.;
    Float_t pzgen=-999.;
    Float_t etagen=-999.;
    Float_t phigen=-999.;
    if(fReadMC){
      TParticle* part = ((AliMCParticle*)mcEvent->GetTrack(TMath::Abs(trlabel)))->Particle();
      ptgen=part->Pt();
      //     pgen=part->P();
      pxgen=part->Px();
      pygen=part->Py();
      pzgen=part->Pz();
      etagen=part->Eta();
      phigen=part->Phi();
    }

    Int_t indexn=0;
    xnt[indexn++]=pttrack;
    xnt[indexn++]=pxtrack;
    xnt[indexn++]=pytrack;
    xnt[indexn++]=pztrack;
    xnt[indexn++]=track->Eta();
    xnt[indexn++]=track->Phi();
    xnt[indexn++]=impactXY;
    xnt[indexn++]=impactZ;
    for(Int_t iLay=0; iLay<4; iLay++) xnt[indexn++]=dedxlay[iLay];
    xnt[indexn++]=dedx;
    xnt[indexn++]=chi2;
    xnt[indexn++]=status;
    xnt[indexn++]=isITSrefit;
    xnt[indexn++]=isTPCin;
    xnt[indexn++]=isTPCrefit;
    xnt[indexn++]=isPureSA;
    xnt[indexn++]=nITSclus;
    xnt[indexn++]=nTPCclus;
    xnt[indexn++]=clumap;
    for(Int_t iLay=0; iLay<6; iLay++) xnt[indexn++]=(Float_t)statusLay[iLay];
    xnt[indexn++]=trlabel;
    xnt[indexn++]=ptgen;
    xnt[indexn++]=pxgen;
    xnt[indexn++]=pygen;
    xnt[indexn++]=pzgen;
    xnt[indexn++]=etagen;
    xnt[indexn++]=phigen;
    xnt[indexn++]=ntracks;
    xnt[indexn++]=ntracklets;
    xnt[indexn++]=centrality;
    xnt[indexn++]=xvert;
    xnt[indexn++]=yvert;
    xnt[indexn++]=zvert;
    xnt[indexn++]=esd->GetEventNumberInFile();

    if(indexn>ntSize) printf("AliAnalysisTaskITSsaTracks: ERROR ntuple insexout of range\n");
    if(fFillNtuple) fNtupleTracks->Fill(xnt);

    if(!(status&AliESDtrack::kITSrefit)) continue;

    Int_t iTrackType=-1;
    if(status&AliESDtrack::kTPCin){
      iTrackType=kTypeTPCITS;
    }else{
      if(status&AliESDtrack::kITSpureSA) iTrackType=kTypeITSpureSA;
      else iTrackType=kTypeITSsa;
    }
    if(iTrackType<0 || iTrackType>=kNtrackTypes) continue;

    Int_t nSPDPoints=0;
    for(Int_t i=0; i<2; i++){
      if(clumap&(1<<i)) ++nSPDPoints;
    }
    Int_t nPointsForPid=0;
    for(Int_t i=2; i<6; i++){
      if(clumap&(1<<i)) ++nPointsForPid;
    }
    if(iTrackType==kTypeTPCITS) nITSTPCtracks+=1;
    if(iTrackType==kTypeITSpureSA) nPureSAtracks+=1;
    Double_t eta=track->Eta();
    Double_t phi=track->Phi();
    fHistEtaPhiAny[iTrackType]->Fill(eta,phi);
    if(nSPDPoints>=1) fHistEtaPhi1SPD[iTrackType]->Fill(eta,phi);
    if(nSPDPoints>=1 && nPointsForPid>=3) fHistEtaPhi4Clu[iTrackType]->Fill(eta,phi);
    if(nITSclus==6) fHistEtaPhi6Clu[iTrackType]->Fill(eta,phi);
    fHistEtaPhiLay[7*iTrackType]->Fill(eta,phi);

    fHistCluInLayVsPt[iTrackType]->Fill(ptrack,-1);
    for(Int_t iBit=0; iBit<6; iBit++){
      if(clumap&(1<<iBit)){
	fHistCluInLayVsPt[iTrackType]->Fill(pttrack,iBit);
	fHistEtaPhiLay[7*iTrackType+iBit+1]->Fill(eta,phi);
      }
    }



    if(nITSclus<fMinITSpts)continue;
    if((chi2/nITSclus) > fMaxITSChi2Clu) continue;
    if(nSPDPoints<fMinSPDpts) continue;

    fHistPt[iTrackType]->Fill(pttrack);
    fHistEtaPhi[iTrackType]->Fill(eta,phi);

    fHistChi2[iTrackType]->Fill(chi2/nITSclus);
    fHistNclu[iTrackType]->Fill(nITSclus);
    if(nPointsForPid==2){
      fHistdedxvsP2cls[iTrackType]->Fill(ptrack,dedx);
    }
    else if(nPointsForPid==3){
      fHistdedxvsP3cls[iTrackType]->Fill(ptrack,dedx);
    }
    else if(nPointsForPid==4){
      fHistdedxvsP4cls[iTrackType]->Fill(ptrack,dedx);
    }
    if(fReadMC){
      if(trlabel<0){
	fHistPtFake[iTrackType]->Fill(pttrack);
	fHistEtaPhiFake[iTrackType]->Fill(eta,phi);
	fHistChi2Fake[iTrackType]->Fill(chi2/nITSclus);
	fHistNcluFake[iTrackType]->Fill(nITSclus);
      }else{
	fHistPtGood[iTrackType]->Fill(pttrack);
	fHistEtaPhiGood[iTrackType]->Fill(eta,phi);
	fHistChi2Good[iTrackType]->Fill(chi2/nITSclus);
	fHistNcluGood[iTrackType]->Fill(nITSclus);
      }
    }

    Int_t hadronSpecie=-1;
    if(fReadMC && fUseMCId){
      if(trlabel>=0){
	TParticle* part = ((AliMCParticle*)mcEvent->GetTrack(trlabel))->Particle();
	Int_t pdg=TMath::Abs(part->GetPdgCode());
	if(pdg==211) hadronSpecie=kPion;
	else if(pdg==321) hadronSpecie=kKaon;
	else if(pdg==2212) hadronSpecie=kProton;
      }
    }else{
      if(nPointsForPid<fMinPtsforPid)continue;
      AliITSPIDResponse pidits(fReadMC);
      Float_t nSigmaPion=pidits.GetNumberOfSigmas(ptrack,dedx,AliPID::kPion,nPointsForPid,isSA);
      if(nSigmaPion>-2. && nSigmaPion<2.){
	hadronSpecie=kPion;
      }else{
	Float_t nSigmaKaon=pidits.GetNumberOfSigmas(ptrack,dedx,AliPID::kKaon,nPointsForPid,isSA);
	if(nSigmaKaon>-2. && nSigmaKaon<2.){
	  hadronSpecie=kKaon;
	}else{
	  Float_t nSigmaProton=pidits.GetNumberOfSigmas(ptrack,dedx,AliPID::kProton,nPointsForPid,isSA);
	  if(nSigmaProton>-2. && nSigmaProton<2.){
	    hadronSpecie=kProton;
	  }
	}
      } 
    }
    if(hadronSpecie<0) continue;
    if(iTrackType==kTypeTPCITS){ // TPC+ITS tracks
      fHistPtTPCITS[hadronSpecie]->Fill(pttrack);
      fHistEtaPhiTPCITS[hadronSpecie]->Fill(eta,phi);
      fHistNcluTPCITS[hadronSpecie]->Fill(pttrack,nITSclus);
      fHistCluInLayTPCITS[hadronSpecie]->Fill(pttrack,-1.);
      for(Int_t iBit=0; iBit<6; iBit++){
	if(clumap&(1<<iBit)) fHistCluInLayTPCITS[hadronSpecie]->Fill(pttrack,iBit);
      }
    }else if(iTrackType==kTypeITSpureSA){ // TPC+ITS tracks
      fHistPtITSpureSA[hadronSpecie]->Fill(pttrack);
      fHistEtaPhiITSpureSA[hadronSpecie]->Fill(eta,phi);
      fHistNcluITSpureSA[hadronSpecie]->Fill(pttrack,nITSclus);
      fHistd0rphiITSpureSA[hadronSpecie]->Fill(pttrack,impactXY);
      fHistd0zITSpureSA[hadronSpecie]->Fill(pttrack,impactZ);
      fHistCluInLayITSpureSA[hadronSpecie]->Fill(pttrack,-1.);
      Int_t outerLay=-1;
      for(Int_t iBit=0; iBit<6; iBit++){
	if(clumap&(1<<iBit)){
	  fHistCluInLayITSpureSA[hadronSpecie]->Fill(pttrack,iBit);
	  if(iBit>outerLay) outerLay=iBit;
	}
	}
      fHistOuterLayITSpureSA[hadronSpecie]->Fill(pttrack,outerLay);	
      
	
      if(fReadMC){  
	if(trlabel>=0){
	  TParticle* part =  ((AliMCParticle*)mcEvent->GetTrack(trlabel))->Particle();
	  ptgen=part->Pt();
	  Float_t invpttrack=track->OneOverPt();
	  Float_t invptgen=0.;
	  if(ptgen>0.) invptgen=1./ptgen;
	  fHistMCPtResid[hadronSpecie]->Fill(pttrack,pttrack-ptgen);
	  fHistMCPtRelResid[hadronSpecie]->Fill(pttrack,(pttrack-ptgen)/ptgen);
	  fHistMCInvPtResid[hadronSpecie]->Fill(pttrack,invpttrack-invptgen);
	  fHistMCInvPtRelResid[hadronSpecie]->Fill(pttrack,(invpttrack-invptgen)/invptgen);	  
	  Float_t deltaphi=track->Phi()-part->Phi();
	  if(deltaphi<-TMath::Pi()) deltaphi+=2*TMath::Pi();
	  if(deltaphi>TMath::Pi()) deltaphi-=2*TMath::Pi();
	  fHistMCPhiResid->Fill(pttrack,deltaphi);
	}
      }
    }else if(iTrackType==kTypeITSsa){ // TPC+ITS tracks
      fHistPtITSsa[hadronSpecie]->Fill(pttrack);
      fHistEtaPhiITSsa[hadronSpecie]->Fill(eta,phi);
      fHistNcluITSsa[hadronSpecie]->Fill(pttrack,nITSclus);
      fHistCluInLayITSsa[hadronSpecie]->Fill(pttrack,-1.);
      for(Int_t iBit=0; iBit<6; iBit++){
	if(clumap&(1<<iBit)) fHistCluInLayITSsa[hadronSpecie]->Fill(pttrack,iBit);
      }
    }
  
    // match TPCITS with ITSpureSA
    if(nITSclus<fMinITSpts || nTPCclus<fMinTPCpts) continue;      
    Int_t idxMI[12],idxSA[12];
    for(Int_t icl=0; icl<12; icl++){ 
      idxMI[icl]=-1; 
      idxSA[icl]=-1;
    }
    Int_t ncls=track->GetClusters(0,idxMI);
    if(fMinITSpts<6){
      Bool_t accept=kTRUE;
      for(Int_t ilay=0; ilay<6; ilay++){
	if(fRequirePoint[ilay]){
	  Int_t mask = 1<<ilay;
	  if(!(clumap & mask)){ 
	    accept=kFALSE;
	    break;
	  }
	}
      }
      if(!accept) continue;
    }
    // Sort
    for(Int_t i=0;i<12;i++){
      for(Int_t j=i+1;j<12;j++){
	if(idxMI[j]>idxMI[i]){
	  Int_t tmp=idxMI[j];
	  idxMI[j]=idxMI[i];
	  idxMI[i]=tmp;
	}
      }
    }
    //    for(Int_t i=0; i<12; i++) printf("%d ",idxMI[i]);
    //    printf("\n");
    if(idxMI[0]<0 && idxMI[0]==idxMI[1]) continue;
    Bool_t matched=kFALSE;
    for (Int_t iTrack2 = 0; iTrack2 < ntracks; iTrack2++) {
      if(matched) break;
      if(iTrack2==iTrack) continue;
      AliESDtrack* track2 = esd->GetTrack(iTrack2);
      Int_t status2=track2->GetStatus();
      if(!(status2&AliESDtrack::kITSrefit)) continue;
      if(!(status2&AliESDtrack::kITSpureSA)) continue;      
      Int_t clumap2=track2->GetITSClusterMap();
      Int_t nITSclus2=track2->GetNcls(0);
      Int_t nTPCclus2=track2->GetNcls(1);
      if(nITSclus2<fMinITSpts || nTPCclus2>0) continue; 
      Int_t ncls2=track2->GetClusters(0,idxSA);
      if(ncls2!=ncls) continue;
      if(fMinITSpts<6){
	Bool_t accept=kTRUE;
	for(Int_t ilay=0; ilay<6; ilay++){
	  if(fRequirePoint[ilay]){
	    Int_t mask = 1<<ilay;
	    if(!(clumap2 & mask)){ 
	      accept=kFALSE;
	      break;
	    }
	  }
	}
	if(!accept) continue;
      }
      // Sort
      for(Int_t i=0;i<12;i++){
	for(Int_t j=i+1;j<12;j++){
	  if(idxSA[j]>idxSA[i]){
	    Int_t tmp=idxSA[j];
	    idxSA[j]=idxSA[i];
	    idxSA[i]=tmp;
	  }
	}
      }
      Int_t match=0;
      for(Int_t icl=0; icl<ncls; icl++){
	if(idxSA[icl]!=idxMI[icl]){
	  match=0; 
	  break;
	}
	else match++;
      }
      if(match==ncls && match>0){
	matched=kTRUE;
	Float_t pt1=track->Pt();
	Float_t pt2=track2->Pt();
	Float_t ptm1=track->OneOverPt();
	Float_t ptm2=track2->OneOverPt();
	fHistPtResid[hadronSpecie]->Fill(pt1,pt2-pt1);
	fHistPtRelResid[hadronSpecie]->Fill(pt1,(pt2-pt1)/pt1);
	fHistInvPtResid[hadronSpecie]->Fill(pt1,ptm2-ptm1);
	fHistInvPtRelResid[hadronSpecie]->Fill(pt1,(ptm2-ptm1)/ptm1);
	Float_t deltaphi=track->Phi()-track2->Phi();
	if(deltaphi<-TMath::Pi()) deltaphi+=2*TMath::Pi();
	if(deltaphi>TMath::Pi()) deltaphi-=2*TMath::Pi();
	fHistPhiResid->Fill(pt1,deltaphi);
      }
    }
  }
  delete esdTrackCutsTPC;
  fHistPureSAtracksVsTracklets->Fill(ntracklets,nPureSAtracks);
  fHistITSTPCtracksVsTracklets->Fill(ntracklets,nITSTPCtracks);
  PostData(1,fOutput);
  
}
//______________________________________________________________________________
void AliAnalysisTaskITSsaTracks::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents= dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  printf("AliAnalysisTaskITSsaTracks::Terminate --- Number of events: read = %.0f  analysed = %.0f\n",fHistNEvents->GetBinContent(1),fHistNEvents->GetBinContent(15));
  return;
}





