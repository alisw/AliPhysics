/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/*
    Comments to be written here: 
    1. What do we calibrate.
    2. How to interpret results
    3. Simple example
    4. Analysis using debug streamers.



    3.Simple example
    // To make cosmic scan the user interaction neccessary
    //
    .x ~/UliStyle.C
    gSystem->Load("libANALYSIS");
    gSystem->Load("libTPCcalib");
    TFile fcalib("CalibObjects.root");
    TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
    AliTPCcalibCosmic * cosmic = ( AliTPCcalibCosmic *)array->FindObject("cosmicTPC");
    

    //
    //
    gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
    gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
    AliXRDPROOFtoolkit tool;
    TChain * chainCosmic = tool.MakeChainRandom("cosmic.txt","Track0",0,10000);
    TChain * chainBudget = tool.MakeChainRandom("cosmic.txt","material",0,1000); 
    TCut cutptV="abs(1/pt0V-1/pt1V)<0.1";
    TCut cutptI="abs(1/pt0In-1/pt1In)<0.5";
    TCut cutncl="nclmin>120";
    TCut cutDz="abs(p0.fP[1])<50";
    TCut cutDr="abs(p0.fP[0])<50";
    //
    chainBudget->Draw(">>listB",cutptV+cutptI+cutncl+cutDr+cutDz,"entryList");
    TEntryList *elistB = (TEntryList*)gDirectory->Get("listB");
    chainBudget->SetEntryList(elistB);
    
    chainBudget->SetAlias("dptrel","(pt0V-pt1V)/((pt0V+pt1V)*0.5)");
    chainBudget->SetAlias("dptInrel","(pt0In-pt1In)/((pt0In+pt1In)*0.5)");
    chainBudget->SetAlias("ptcorr","(pt0In-pt0V)/(pt0V)+(pt1V-pt1In)/(pt1In)");
*/



#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TDatabasePDG.h"

#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"

#include "AliTracker.h"
#include "AliMagF.h"
#include "AliTPCCalROC.h"

#include "AliLog.h"

#include "AliTPCcalibCosmic.h"
#include "TTreeStream.h"
#include "AliTPCTracklet.h"
//#include "AliESDcosmic.h"


ClassImp(AliTPCcalibCosmic)


AliTPCcalibCosmic::AliTPCcalibCosmic() 
  :AliTPCcalibBase(),
   fHistNTracks(0),
   fClusters(0),
   fModules(0),
   fHistPt(0),
   fDeDx(0),
   fDeDxMIP(0),
   fMIPvalue(1), 
   fCutMaxD(5),        // maximal distance in rfi ditection
   fCutMaxDz(40),      // maximal distance in z ditection
   fCutTheta(0.03),    // maximal distan theta
   fCutMinDir(-0.99)   // direction vector products
{  
  AliInfo("Default Constructor");    
  for (Int_t ihis=0; ihis<6;ihis++){
    fHistoDelta[ihis]=0;
    fHistoPull[ihis]=0;
  }
  for (Int_t ihis=0; ihis<4;ihis++){
    fHistodEdxMax[ihis]    =0;
    fHistodEdxTot[ihis]    =0;
  }
}


AliTPCcalibCosmic::AliTPCcalibCosmic(const Text_t *name, const Text_t *title) 
  :AliTPCcalibBase(),
   fHistNTracks(0),
   fClusters(0),
   fModules(0),
   fHistPt(0),
   fDeDx(0),
   fDeDxMIP(0),
   fMIPvalue(1),
   fCutMaxD(5),        // maximal distance in rfi ditection 
   fCutMaxDz(40),      // maximal distance in z ditection
   fCutTheta(0.03),    // maximal distan theta
   fCutMinDir(-0.99)   // direction vector products
{  
  SetName(name);
  SetTitle(title);

  fHistNTracks = new TH1F("ntracks","Number of Tracks per Event; number of tracks per event; number of tracks",501,-0.5,500.5);
  fClusters = new TH1F("signal","Number of Clusters per track; number of clusters per track n_{cl}; counts",160,0,160);
  fModules = new TH2F("sector","Acorde hits; z (cm); x(cm)",1200,-650,650,600,-700,700);
  fHistPt = new TH1F("Pt","Pt distribution; p_{T} (GeV); counts",2000,0,50);
  fDeDx = new TH2F("DeDx","dEdx; momentum p (GeV); TPC signal (a.u.)",500,0.01,100.,500,2.,1000);
  BinLogX(fDeDx);
  fDeDxMIP =  new TH1F("DeDxMIP","MIP region; TPC signal (a.u.);counts ",500,2.,1000);
  Init();
  AliInfo("Non Default Constructor");  
  //
}

AliTPCcalibCosmic::~AliTPCcalibCosmic(){
  //
  //
  //
  for (Int_t ihis=0; ihis<6;ihis++){
    delete fHistoDelta[ihis];
    delete fHistoPull[ihis];
  }
  for (Int_t ihis=0; ihis<4;ihis++){
    delete fHistodEdxTot[ihis];
    delete fHistodEdxMax[ihis];
  }

  delete fHistNTracks;            //  histogram showing number of ESD tracks per event
  delete fClusters;               //  histogram showing the number of clusters per track
  delete fModules;                //  2d histogram of tracks which are propagated to the ACORDE scintillator array
  delete fHistPt;                 //  Pt histogram of reconstructed tracks
  delete fDeDx;                   //  dEdx spectrum showing the different particle types
  delete fDeDxMIP;                //  TPC signal close to the MIP region of muons 0.4 < p < 0.45 GeV
}


void AliTPCcalibCosmic::Init(){
  //
  // init component
  // Make performance histograms
  //

  // tracking performance bins
  // 0 - delta of interest
  // 1 - min (track0, track1) number of clusters
  // 2 - R  - vertex radius
  // 3 - P1 - mean z
  // 4 - P2 - snp(phi)    at inner wall of TPC
  // 5 - P3 - tan(theta)  at inner wall of TPC
  // 6 - P4 - 1/pt mean
  // 7 - pt - pt mean
  // 8 - alpha

  Double_t xminTrack[9], xmaxTrack[9];
  Int_t binsTrack[9];
  TString axisName[9];
  //
  binsTrack[0] =100;
  axisName[0]  ="#Delta";
  //
  binsTrack[1] =8;
  xminTrack[1] =80; xmaxTrack[1]=160;
  axisName[1]  ="N_{cl}";
  //
  binsTrack[2] =10;
  xminTrack[2] =0; xmaxTrack[2]=90;  // 
  axisName[2]  ="dca_{r} (cm)";
  //
  binsTrack[3] =25;
  xminTrack[3] =-250; xmaxTrack[3]=250;  // 
  axisName[3]  ="z (cm)";
  //
  binsTrack[4] =10;
  xminTrack[4] =-0.8; xmaxTrack[4]=0.8;  // 
  axisName[4]  ="sin(#phi)";
  //
  binsTrack[5] =10;
  xminTrack[5] =-1; xmaxTrack[5]=1;  // 
  axisName[5]  ="tan(#theta)";
  //
  binsTrack[6] =40;
  xminTrack[6] =-2; xmaxTrack[6]=2;  // 
  axisName[6]  ="1/pt (1/GeV)";
  //
  binsTrack[7] =40;
  xminTrack[7] =0.2; xmaxTrack[7]=50;  // 
  axisName[7]  ="pt (GeV)";
  //
  binsTrack[8] =18;
  xminTrack[8] =0; xmaxTrack[8]=TMath::Pi();  // 
  axisName[8]  ="alpha";
  //
  // delta y
  xminTrack[0] =-1; xmaxTrack[0]=1;  // 
  fHistoDelta[0] = new THnSparseS("#Delta_{Y} (cm)","#Delta_{Y} (cm)", 9, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fHistoPull[0] = new THnSparseS("#Delta_{Y} (unit)","#Delta_{Y} (unit)", 9, binsTrack,xminTrack, xmaxTrack);
  //
  // delta z
  xminTrack[0] =-1; xmaxTrack[0]=1;  // 
  fHistoDelta[1] = new THnSparseS("#Delta_{Z} (cm)","#Delta_{Z} (cm)", 9, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fHistoPull[1] = new THnSparseS("#Delta_{Z} (unit)","#Delta_{Z} (unit)", 9, binsTrack,xminTrack, xmaxTrack);
  //
  // delta P2
  xminTrack[0] =-10; xmaxTrack[0]=10;  // 
  fHistoDelta[2] = new THnSparseS("#Delta_{#phi} (mrad)","#Delta_{#phi} (mrad)", 9, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fHistoPull[2] = new THnSparseS("#Delta_{#phi} (unit)","#Delta_{#phi} (unit)", 9, binsTrack,xminTrack, xmaxTrack);
  //
  // delta P3
  xminTrack[0] =-10; xmaxTrack[0]=10;  // 
  fHistoDelta[3] = new THnSparseS("#Delta_{#theta} (mrad)","#Delta_{#theta} (mrad)", 9, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fHistoPull[3] = new THnSparseS("#Delta_{#theta} (unit)","#Delta_{#theta} (unit)", 9, binsTrack,xminTrack, xmaxTrack);
  //
  // delta P4
  xminTrack[0] =-0.2; xmaxTrack[0]=0.2;  // 
  fHistoDelta[4] = new THnSparseS("#Delta_{1/pt} (1/GeV)","#Delta_{1/pt} (1/GeV)", 9, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fHistoPull[4] = new THnSparseS("#Delta_{1/pt} (unit)","#Delta_{1/pt} (unit)", 9, binsTrack,xminTrack, xmaxTrack);
  
  //
  // delta Pt
  xminTrack[0] =-0.5; xmaxTrack[0]=0.5;  // 
  fHistoDelta[5] = new THnSparseS("#Delta_{pt}/p_{t}","#Delta_{pt}/p_{t}", 9, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fHistoPull[5] = new THnSparseS("#Delta_{pt}/p_{t} (unit)","#Delta_{pt}/p_{t} (unit)", 9, binsTrack,xminTrack, xmaxTrack);
  //

  for (Int_t idedx=0;idedx<4;idedx++){
    xminTrack[0] =0.5; xmaxTrack[0]=1.5;  // 
    binsTrack[1] =40;
    xminTrack[1] =10; xmaxTrack[1]=160;

    fHistodEdxMax[idedx] = new THnSparseS(Form("dEdx_{MaxUp}/dEdx_{MaxDown}_Pad%d",idedx),
					  Form("dEdx_{MaxUp}/dEdx_{MaxDown}_Pad%d",idedx), 
					  9, binsTrack,xminTrack, xmaxTrack);
    fHistodEdxTot[idedx] = new THnSparseS(Form("dEdx_{TotUp}/dEdx_{TotDown}_Pad%d",idedx),
					  Form("dEdx_{TotUp}/dEdx_{TotDown}_Pad%d",idedx), 
					  9, binsTrack,xminTrack, xmaxTrack);
  }
  


  for (Int_t ivar=0;ivar<6;ivar++){
    for (Int_t ivar2=0;ivar2<9;ivar2++){      
      fHistoDelta[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
      fHistoDelta[ivar]->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
      fHistoPull[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
      fHistoPull[ivar]->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
      BinLogX(fHistoDelta[ivar],7);
      BinLogX(fHistoPull[ivar],7);
      if (ivar<4){
	fHistodEdxMax[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
	fHistodEdxMax[ivar]->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
	fHistodEdxTot[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
	fHistodEdxTot[ivar]->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
	BinLogX(fHistodEdxMax[ivar],7);
	BinLogX(fHistodEdxTot[ivar],7);
      }
    }
  }
}


void AliTPCcalibCosmic::Add(const AliTPCcalibCosmic* cosmic){
  //
  //
  //
  for (Int_t ivar=0; ivar<6;ivar++){
    if (fHistoDelta[ivar] && cosmic->fHistoDelta[ivar]){
      fHistoDelta[ivar]->Add(cosmic->fHistoDelta[ivar]);
    }
    if (fHistoPull[ivar] && cosmic->fHistoPull[ivar]){
      fHistoPull[ivar]->Add(cosmic->fHistoPull[ivar]);
    }
  }
  for (Int_t ivar=0; ivar<4;ivar++){
    if (fHistodEdxMax[ivar] && cosmic->fHistodEdxMax[ivar]){
      fHistodEdxMax[ivar]->Add(cosmic->fHistodEdxMax[ivar]);
    }
    if (fHistodEdxTot[ivar] && cosmic->fHistodEdxTot[ivar]){
      fHistodEdxTot[ivar]->Add(cosmic->fHistodEdxTot[ivar]);
    }
  }
}




void AliTPCcalibCosmic::Process(AliESDEvent *event) {
  //
  //
  //
  if (!event) {
    Printf("ERROR: ESD not available");
    return;
  }  
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!esdFriend) {
   Printf("ERROR: esdFriend not available");
   return;
  }
   

  FindPairs(event); // nearly everything takes place in find pairs...

  if (GetDebugLevel()>20) printf("Hallo world: Im here and processing an event\n");
  Int_t ntracks=event->GetNumberOfTracks(); 
  fHistNTracks->Fill(ntracks);
  if (ntracks==0) return;
 //  AliESDcosmic cosmicESD;    
//   TTreeSRedirector * cstream =  GetDebugStreamer();
//   cosmicESD.SetDebugStreamer(cstream);
//   cosmicESD.ProcessEvent(event);
//   if (cstream) cosmicESD.DumpToTree();
      
  
}


void AliTPCcalibCosmic::FillHistoPerformance(const AliExternalTrackParam *par0, const AliExternalTrackParam *par1, const AliExternalTrackParam *inner0, const AliExternalTrackParam */*inner1*/, AliTPCseed *seed0,  AliTPCseed *seed1, const AliExternalTrackParam *param0Combined ){
  //
  // par0,par1       - parameter of tracks at DCA 0
  // inner0,inner1   - parameter of tracks at the TPC entrance
  // seed0, seed1    - detailed track information
  // param0Combined  - Use combined track parameters for binning
  //
  Int_t kMinCldEdx =20;
  Int_t ncl0 = seed0->GetNumberOfClusters();
  Int_t ncl1 = seed1->GetNumberOfClusters();

  const Double_t kpullCut    = 10;
  Double_t x[9];
  Double_t xyz0[3];
  Double_t xyz1[3];
  par0->GetXYZ(xyz0);
  par1->GetXYZ(xyz1);
  Double_t radius0 = TMath::Sqrt(xyz0[0]*xyz0[0]+xyz0[1]*xyz0[1]);
  Double_t radius1 = TMath::Sqrt(xyz1[0]*xyz1[0]+xyz1[1]*xyz1[1]);
  inner0->GetXYZ(xyz0);
  Double_t alpha = TMath::ATan2(xyz0[1],xyz0[0]);
  // bin parameters
  x[1] = TMath::Min(ncl0,ncl1);
  x[2] = (radius0+radius1)*0.5;
  x[3] = param0Combined->GetZ();
  x[4] = inner0->GetSnp();
  x[5] = param0Combined->GetTgl();
  x[6] = param0Combined->GetSigned1Pt();
  x[7] = param0Combined->Pt();
  x[8] = alpha;
  // deltas
  Double_t delta[6];
  Double_t sigma[6];
  delta[0] = (par0->GetY()+par1->GetY());
  delta[1] = (par0->GetZ()-par1->GetZ());
  delta[2] = (par0->GetAlpha()-par1->GetAlpha()-TMath::Pi());
  delta[3] = (par0->GetTgl()+par1->GetTgl());
  delta[4] = (par0->GetParameter()[4]+par1->GetParameter()[4]);
  delta[5] = (par0->Pt()-par1->Pt())/((par0->Pt()+par1->Pt())*0.5);
  //
  sigma[0] = TMath::Sqrt(par0->GetSigmaY2()+par1->GetSigmaY2());
  sigma[1] = TMath::Sqrt(par0->GetSigmaZ2()+par1->GetSigmaZ2());
  sigma[2] = TMath::Sqrt(par0->GetSigmaSnp2()+par1->GetSigmaSnp2());
  sigma[3] = TMath::Sqrt(par0->GetSigmaTgl2()+par1->GetSigmaTgl2());
  sigma[4] = TMath::Sqrt(par0->GetSigma1Pt2()+par1->GetSigma1Pt2());
  sigma[5] = sigma[4]*((par0->Pt()+par1->Pt())*0.5);
  //
  Bool_t isOK = kTRUE;
  for (Int_t ivar=0;ivar<6;ivar++){
    if (sigma[ivar]==0) isOK=kFALSE;
    x[0]= delta[ivar]/sigma[ivar];
    if (TMath::Abs(x[0])>kpullCut) isOK = kFALSE;
  }
  //

  if (isOK) for (Int_t ivar=0;ivar<6;ivar++){
    x[0]= delta[ivar]/TMath::Sqrt(2);
    if (ivar==2 || ivar ==3) x[0]*=1000;
    fHistoDelta[ivar]->Fill(x);
    if (sigma[ivar]>0){
      x[0]= delta[ivar]/sigma[ivar];
      fHistoPull[ivar]->Fill(x);
    }
  }

  //						
  // Fill dedx performance
  //
  for (Int_t ipad=0; ipad<4;ipad++){
    //
    //
    //
    Int_t row0=0;
    Int_t row1=160;
    if (ipad==0) row1=63;
    if (ipad==1) {row0=63; row1=63+64;}
    if (ipad==2) {row0=128;}
    Int_t   nclUp       = TMath::Nint(seed0->CookdEdxAnalytical(0.01,0.7,0,row0,row1,2));
    Int_t   nclDown     = TMath::Nint(seed1->CookdEdxAnalytical(0.01,0.7,0,row0,row1,2));
    Int_t   minCl       = TMath::Min(nclUp,nclDown);
    if (minCl<kMinCldEdx) continue;
    x[1] = minCl;
    //
    Float_t dEdxTotUp   = seed0->CookdEdxAnalytical(0.01,0.7,0,row0,row1);
    Float_t dEdxTotDown = seed1->CookdEdxAnalytical(0.01,0.7,0,row0,row1);
    Float_t dEdxMaxUp   = seed0->CookdEdxAnalytical(0.01,0.7,1,row0,row1);
    Float_t dEdxMaxDown = seed1->CookdEdxAnalytical(0.01,0.7,1,row0,row1);
    //
    if (dEdxTotDown<=0) continue;
    if (dEdxMaxDown<=0) continue;
    x[0]=dEdxTotUp/dEdxTotDown;
    fHistodEdxTot[ipad]->Fill(x);
    x[0]=dEdxMaxUp/dEdxMaxDown;
    fHistodEdxMax[ipad]->Fill(x);
  }


  
}

void AliTPCcalibCosmic::MaterialBudgetDump(AliExternalTrackParam *const par0, AliExternalTrackParam *const par1, const AliExternalTrackParam *inner0, const AliExternalTrackParam *inner1, AliTPCseed *const seed0,  AliTPCseed *const seed1, AliExternalTrackParam *const param0Combined, AliExternalTrackParam *const param1Combined){
  //
  // matrial budget AOD dump
  //
  // par0,par1       - parameter of tracks at DCA 0
  // inner0,inner1   - parameter of tracks at the TPC entrance
  // seed0, seed1    - detailed track information
  // param0Combined  - Use combined track parameters for binning
  // param1Combined  - 
  Double_t p0In = inner0->GetP(); 
  Double_t p1In = inner1->GetP(); 
  Double_t p0V  = par0->GetP(); 
  Double_t p1V  = par1->GetP(); 
  //
  Double_t pt0In = inner0->Pt(); 
  Double_t pt1In = inner1->Pt(); 
  Double_t pt0V  = par0->Pt(); 
  Double_t pt1V  = par1->Pt(); 
  Int_t ncl0 = seed0->GetNumberOfClusters();
  Int_t ncl1 = seed1->GetNumberOfClusters();
  Int_t nclmin=TMath::Min(ncl0,ncl1);
  Double_t sign = (param0Combined->GetSigned1Pt()>0) ? 1:-1.;
  //
  TTreeSRedirector * pcstream =  GetDebugStreamer();
  if (pcstream){
    (*pcstream)<<"material"<<
      "run="<<fRun<<              //  run number
      "event="<<fEvent<<          //  event number
      "time="<<fTime<<            //  time stamp of event
      "trigger="<<fTrigger<<      //  trigger
      "triggerClass="<<&fTriggerClass<<      //  trigger
      "mag="<<fMagF<<             //  magnetic field
      "sign="<<sign<<             // sign of the track
      //
      "ncl0="<<ncl0<<
      "ncl1="<<ncl1<<
      "nclmin="<<nclmin<<
      //
      "p0In="<<p0In<<
      "p1In="<<p1In<<
      "p0V="<<p0V<<
      "p1V="<<p1V<<
      "pt0In="<<pt0In<<
      "pt1In="<<pt1In<<
      "pt0V="<<pt0V<<
      "pt1V="<<pt1V<<
      "p0.="<<par0<<
      "p1.="<<par1<<
      "up0.="<<param0Combined<<
      "up1.="<<param1Combined<<
      "\n";
  }
  
}

void AliTPCcalibCosmic::Analyze() {

  fMIPvalue = CalculateMIPvalue(fDeDxMIP);

  return;

}



void AliTPCcalibCosmic::FindPairs(AliESDEvent *event) {
  //
  // Find cosmic pairs
  // 
  // Track0 is choosen in upper TPC part
  // Track1 is choosen in lower TPC part
  //
  if (GetDebugLevel()>20) printf("Hallo world: Im here\n");
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  Int_t ntracks=event->GetNumberOfTracks(); 
  TObjArray  tpcSeeds(ntracks);
  if (ntracks==0) return;
  Double_t vtxx[3]={0,0,0};
  Double_t svtxx[3]={0.000001,0.000001,100.};
  AliESDVertex vtx(vtxx,svtxx);
  //
  //track loop
  //
  for (Int_t i=0;i<ntracks;++i) {
   AliESDtrack *track = event->GetTrack(i);
   fClusters->Fill(track->GetTPCNcls()); 
  
   const AliExternalTrackParam * trackIn = track->GetInnerParam();
   const AliExternalTrackParam * trackOut = track->GetOuterParam();
   if (!trackIn) continue;
   if (!trackOut) continue;
   if (ntracks>4 && TMath::Abs(trackIn->GetTgl())<0.0015) continue;  // filter laser 


   AliESDfriendTrack *friendTrack = esdFriend->GetTrack(i);
   TObject *calibObject;
   AliTPCseed *seed = 0;
   for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
     if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
   }
   if (seed) tpcSeeds.AddAt(seed,i);

   Double_t meanP = 0.5*(trackIn->GetP() + trackOut->GetP());
   if (seed && track->GetTPCNcls() > 80 + 60/(1+TMath::Exp(-meanP+5))) {
     fDeDx->Fill(meanP, seed->CookdEdxNorm(0.0,0.45,0,0,159));
     //
     if (meanP > 0.4 && meanP < 0.45) fDeDxMIP->Fill(seed->CookdEdxNorm(0.0,0.45,0,0,159));
     //
     if (GetDebugLevel()>0&&meanP>0.2&&seed->CookdEdxNorm(0.0,0.45,0,0,159)>300) {
       TFile *curfile = AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile();
       //if (curfile) printf(">>> p+ in file: %s \t event: %i \t Number of ESD tracks: %i \n", curfile->GetName(), (int)event->GetEventNumberInFile(), (int)ntracks);
       // if (track->GetOuterParam()->GetAlpha()<0) cout << " Polartiy: " << track->GetSign() << endl;
     }

   }

  }

  if (ntracks<2) return;
  //
  // Find pairs
  //
  for (Int_t i=0;i<ntracks;++i) {
    AliESDtrack *track0 = event->GetTrack(i);     
    // track0 - choosen upper part
    if (!track0) continue;
    if (!track0->GetOuterParam()) continue;
    if (track0->GetOuterParam()->GetAlpha()<0) continue;
    Double_t dir0[3];
    track0->GetDirection(dir0);    
    for (Int_t j=0;j<ntracks;++j) {
      if (i==j) continue;
      AliESDtrack *track1 = event->GetTrack(j);   
      //track 1 lower part
      if (!track1) continue;
      if (!track1->GetOuterParam()) continue;
      if (track1->GetOuterParam()->GetAlpha()>0) continue;
      //
      Double_t dir1[3];
      track1->GetDirection(dir1);
      
      AliTPCseed * seed0 = (AliTPCseed*) tpcSeeds.At(i);
      AliTPCseed * seed1 = (AliTPCseed*) tpcSeeds.At(j);
      if (! seed0) continue;
      if (! seed1) continue;
      Float_t dedx0 = seed0->CookdEdxNorm(0.05,0.55,0,0,159);
      Float_t dedx1 = seed1->CookdEdxNorm(0.05,0.55,0,0,159);
      //
      Float_t dedx0I = seed0->CookdEdxNorm(0.05,0.55,0,0,63);
      Float_t dedx1I = seed1->CookdEdxNorm(0.05,0.55,0,0,63);
      //
      Float_t dedx0O = seed0->CookdEdxNorm(0.05,0.55,0,64,159);
      Float_t dedx1O = seed1->CookdEdxNorm(0.05,0.55,0,64,159);
      //
      Float_t dir = (dir0[0]*dir1[0] + dir0[1]*dir1[1] + dir0[2]*dir1[2]);
      Float_t d0  = track0->GetLinearD(0,0);
      Float_t d1  = track1->GetLinearD(0,0);
      //
      // conservative cuts - convergence to be guarantied
      // applying before track propagation
      if (TMath::Abs(d0+d1)>fCutMaxD) continue;   // distance to the 0,0
      if (dir>fCutMinDir) continue;               // direction vector product
      Float_t bz = AliTracker::GetBz();
      Float_t dvertex0[2];   //distance to 0,0
      Float_t dvertex1[2];   //distance to 0,0 
      track0->GetDZ(0,0,0,bz,dvertex0);
      track1->GetDZ(0,0,0,bz,dvertex1);
      if (TMath::Abs(dvertex0[1])>250) continue;
      if (TMath::Abs(dvertex1[1])>250) continue;
      //
      //
      //
      Float_t dmax = TMath::Max(TMath::Abs(d0),TMath::Abs(d1));
      AliExternalTrackParam param0(*track0);
      AliExternalTrackParam param1(*track1);
      //
      // Propagate using Magnetic field and correct fo material budget
      //
      Double_t sign0=-1;
      Double_t sign1=1;
      Double_t maxsnp=0.90;
      AliTracker::PropagateTrackToBxByBz(&param0,dmax+1,TDatabasePDG::Instance()->GetParticle("e-")->Mass(),3,kTRUE,maxsnp,sign0);
      AliTracker::PropagateTrackToBxByBz(&param1,dmax+1,TDatabasePDG::Instance()->GetParticle("e-")->Mass(),3,kTRUE,maxsnp,sign1);
      //
      // Propagate rest to the 0,0 DCA - z should be ignored
      //
      Bool_t b0 = param0.PropagateToDCA(&vtx,bz,1000);
      Bool_t b1 = param1.PropagateToDCA(&vtx,bz,1000);
      //      
      param0.GetDZ(0,0,0,bz,dvertex0);
      param1.GetDZ(0,0,0,bz,dvertex1);
      if (TMath::Abs(param0.GetZ()-param1.GetZ())>fCutMaxDz) continue;
      //
      Double_t xyz0[3];//,pxyz0[3];
      Double_t xyz1[3];//,pxyz1[3];
      param0.GetXYZ(xyz0);
      param1.GetXYZ(xyz1);
      Bool_t isPair = IsPair(&param0,&param1);
      //
      if (isPair) FillAcordeHist(track0);
      //
      // combined track params 
      //
      AliExternalTrackParam *par0U=MakeCombinedTrack(&param0,&param1);
      AliExternalTrackParam *par1U=MakeCombinedTrack(&param1,&param0);


      //
      if (fStreamLevel>0){
	TTreeSRedirector * cstream =  GetDebugStreamer();
	//printf("My stream=%p\n",(void*)cstream);
	AliExternalTrackParam *ip0 = (AliExternalTrackParam *)track0->GetInnerParam();
	AliExternalTrackParam *ip1 = (AliExternalTrackParam *)track1->GetInnerParam();
	AliExternalTrackParam *op0 = (AliExternalTrackParam *)track0->GetOuterParam();
	AliExternalTrackParam *op1 = (AliExternalTrackParam *)track1->GetOuterParam();
	Bool_t isCrossI = ip0->GetZ()*ip1->GetZ()<0;
	Bool_t isCrossO = op0->GetZ()*op1->GetZ()<0;
	Double_t alpha0 = TMath::ATan2(dir0[1],dir0[0]);
	Double_t alpha1 = TMath::ATan2(dir1[1],dir1[0]);
	//
	//
	//
	FillHistoPerformance(&param0, &param1, ip0, ip1, seed0, seed1,par0U);
	MaterialBudgetDump(&param0, &param1, ip0, ip1, seed0, seed1,par0U,par1U);
	if (cstream) {
	  (*cstream) << "Track0" <<
	    "run="<<fRun<<              //  run number
	    "event="<<fEvent<<          //  event number
	    "time="<<fTime<<            //  time stamp of event
	    "trigger="<<fTrigger<<      //  trigger
	    "triggerClass="<<&fTriggerClass<<      //  trigger
	    "mag="<<fMagF<<             //  magnetic field
	    "dir="<<dir<<               //  direction
	    "OK="<<isPair<<             //  will be accepted
	    "b0="<<b0<<                 //  propagate status
	    "b1="<<b1<<                 //  propagate status
	    "crossI="<<isCrossI<<       //  cross inner
	    "crossO="<<isCrossO<<       //  cross outer
	    //
	    "Orig0.=" << track0 <<      //  original track  0
	    "Orig1.=" << track1 <<      //  original track  1
	    "Tr0.="<<&param0<<          //  track propagated to the DCA 0,0
	    "Tr1.="<<&param1<<          //  track propagated to the DCA 0,0	   
	    "Ip0.="<<ip0<<              //  inner param - upper
	    "Ip1.="<<ip1<<              //  inner param - lower
	    "Op0.="<<op0<<              //  outer param - upper
	    "Op1.="<<op1<<              //  outer param - lower
	    "Up0.="<<par0U<<           //  combined track 0
	    "Up1.="<<par1U<<           //  combined track 1
	    //
	    "v00="<<dvertex0[0]<<       //  distance using kalman
	    "v01="<<dvertex0[1]<<       // 
	    "v10="<<dvertex1[0]<<       //
	    "v11="<<dvertex1[1]<<       // 
	    "d0="<<d0<<                 //  linear distance to 0,0
	    "d1="<<d1<<                 //  linear distance to 0,0
	    //
	    //
	    //
	    "x00="<<xyz0[0]<<           // global position close to vertex
	    "x01="<<xyz0[1]<<
	    "x02="<<xyz0[2]<<
	    //
	    "x10="<<xyz1[0]<<           // global position close to vertex
	    "x11="<<xyz1[1]<<
	    "x12="<<xyz1[2]<<
	    //
	    "alpha0="<<alpha0<<
	    "alpha1="<<alpha1<<
	    "dir00="<<dir0[0]<<           // direction upper
	    "dir01="<<dir0[1]<<
	    "dir02="<<dir0[2]<<
	    //
	    "dir10="<<dir1[0]<<           // direction lower
	    "dir11="<<dir1[1]<<
	    "dir12="<<dir1[2]<<
	    //
	    //
	    "Seed0.=" << seed0 <<       //  original seed 0
	    "Seed1.=" << seed1 <<       //  original seed 1
	    //
	    "dedx0="<<dedx0<<           //  dedx0 - all
	    "dedx1="<<dedx1<<           //  dedx1 - all
	    //
	    "dedx0I="<<dedx0I<<         //  dedx0 - inner ROC
	    "dedx1I="<<dedx1I<<         //  dedx1 - inner ROC
	    //
	    "dedx0O="<<dedx0O<<         //  dedx0 - outer ROC
	    "dedx1O="<<dedx1O<<         //  dedx1 - outer ROC
	    "\n";
	}
      }
      delete par0U;
      delete par1U;
    }
  }  
}    




void  AliTPCcalibCosmic::FillAcordeHist(AliESDtrack *upperTrack) {

  // Pt cut to select straight tracks which can be easily propagated to ACORDE which is outside the magnetic field
  if (upperTrack->Pt() < 10 || upperTrack->GetTPCNcls() < 80) return;
    
  const Double_t acordePlane = 850.; // distance of the central Acorde detectors to the beam line at y =0
  const Double_t roof = 210.5;       // distance from x =0 to end of magnet roof

  Double_t r[3];
  upperTrack->GetXYZ(r);
  Double_t d[3];
  upperTrack->GetDirection(d);
  Double_t x,z;
  z = r[2] + (d[2]/d[1])*(acordePlane - r[1]);
  x = r[0] + (d[0]/d[1])*(acordePlane - r[1]);
  
  if (x > roof) {
    x = r[0] + (d[0]/(d[0]+d[1]))*(acordePlane+roof-r[0]-r[1]);
    z = r[2] + (d[2]/(d[0]+d[1]))*(acordePlane+roof-r[0]-r[1]);
  }
  if (x < -roof) {
    x = r[0] + (d[0]/(d[1]-d[0]))*(acordePlane+roof+r[0]-r[1]);	      
    z = r[2] + (d[2]/(d[1]-d[0]))*(acordePlane+roof+r[0]-r[1]);
  } 

  fModules->Fill(z, x);
 
}



Long64_t AliTPCcalibCosmic::Merge(TCollection *const li) {

  TIterator* iter = li->MakeIterator();
  AliTPCcalibCosmic* cal = 0;

  while ((cal = (AliTPCcalibCosmic*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibCosmic::Class())) {
      //Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }
    
    fHistNTracks->Add(cal->GetHistNTracks());
    fClusters->Add(cal-> GetHistClusters());
    fModules->Add(cal->GetHistAcorde());
    fHistPt->Add(cal->GetHistPt());
    fDeDx->Add(cal->GetHistDeDx());
    fDeDxMIP->Add(cal->GetHistMIP());
    Add(cal);
  }
  return 0;
  
}


Bool_t  AliTPCcalibCosmic::IsPair(AliExternalTrackParam *tr0, AliExternalTrackParam *tr1){
  //
  //
  /*
  // 0. Same direction - OPOSITE  - cutDir +cutT    
  TCut cutDir("cutDir","dir<-0.99")
  // 1. 
  TCut cutT("cutT","abs(Tr1.fP[3]+Tr0.fP[3])<0.03")
  //
  // 2. The same rphi 
  TCut cutD("cutD","abs(Tr0.fP[0]+Tr1.fP[0])<5")
  //
  //
  //
  TCut cutPt("cutPt","abs(Tr1.fP[4]+Tr0.fP[4])<1&&abs(Tr0.fP[4])+abs(Tr1.fP[4])<10");  
  // 1/Pt diff cut
  */
  const Double_t *p0 = tr0->GetParameter();
  const Double_t *p1 = tr1->GetParameter();
  if (TMath::Abs(p0[3]+p1[3])>fCutTheta) return kFALSE;
  if (TMath::Abs(p0[1]-p1[1])>fCutMaxDz) return kFALSE;
  if (TMath::Abs(p0[0]+p1[0])>fCutMaxD)  return kFALSE;
  
  Double_t d0[3], d1[3];
  tr0->GetDirection(d0);    
  tr1->GetDirection(d1);       
  if (d0[0]*d1[0] + d0[1]*d1[1] + d0[2]*d1[2] >fCutMinDir) return kFALSE;
  //
  return kTRUE;  
}
 


Double_t AliTPCcalibCosmic::CalculateMIPvalue(TH1F * hist) {

  TF1 * funcDoubleGaus = new TF1("funcDoubleGaus", "gaus(0)+gaus(3)",0,1000);
  funcDoubleGaus->SetParameters(hist->GetEntries()*0.75,hist->GetMean()/1.3,hist->GetMean()*0.10,
				hist->GetEntries()*0.25,hist->GetMean()*1.3,hist->GetMean()*0.10);
  hist->Fit(funcDoubleGaus);
  Double_t aMIPvalue = TMath::Min(funcDoubleGaus->GetParameter(1),funcDoubleGaus->GetParameter(4));

  delete funcDoubleGaus;

  return aMIPvalue;

}




void AliTPCcalibCosmic::CalculateBetheParams(TH2F */*hist*/, Double_t * /*initialParam*/) {
  //
  // Not implemented yet
  //
  return;

}


void AliTPCcalibCosmic::BinLogX(THnSparse *const h, Int_t axisDim) {

  // Method for the correct logarithmic binning of histograms

  TAxis *axis = h->GetAxis(axisDim);
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];

  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);

  for (int i = 1; i <= bins; i++) {
   newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete newBins;

}


void AliTPCcalibCosmic::BinLogX(TH1 *const h) {

  // Method for the correct logarithmic binning of histograms

  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
   
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete newBins;
  
}


AliExternalTrackParam *AliTPCcalibCosmic::MakeTrack(const AliExternalTrackParam *track0, const AliExternalTrackParam *track1){
  //
  // 
  //
  AliExternalTrackParam *par1R= new AliExternalTrackParam(*track1);
  par1R->Rotate(track0->GetAlpha());
  par1R->PropagateTo(track0->GetX(),AliTracker::GetBz()); 
  //
  //
  Double_t * param = (Double_t*)par1R->GetParameter();
  Double_t * covar = (Double_t*)par1R->GetCovariance();

  param[0]*=1;  //OK
  param[1]*=1;  //OK
  param[2]*=1;  //?
  param[3]*=-1; //OK
  param[4]*=-1; //OK
  //
  covar[6] *=-1.; covar[7] *=-1.; covar[8] *=-1.;
  //covar[10]*=-1.; covar[11]*=-1.; covar[12]*=-1.;
  covar[13]*=-1.;
  return par1R;
}

AliExternalTrackParam *AliTPCcalibCosmic::MakeCombinedTrack(const AliExternalTrackParam *track0, const AliExternalTrackParam *track1){
  //
  // Make combined track
  //
  //
  AliExternalTrackParam * par1T = MakeTrack(track0,track1);
  AliExternalTrackParam * par0U = new AliExternalTrackParam(*track0);
  //
  UpdateTrack(*par0U,*par1T);
  delete par1T;
  return par0U;
}


void AliTPCcalibCosmic::UpdateTrack(AliExternalTrackParam &track1, const AliExternalTrackParam &track2){
  //
  // Update track 1 with track 2
  //
  //
  //
  TMatrixD vecXk(5,1);    // X vector
  TMatrixD covXk(5,5);    // X covariance 
  TMatrixD matHk(5,5);    // vector to mesurement
  TMatrixD measR(5,5);    // measurement error 
  TMatrixD vecZk(5,1);    // measurement
  //
  TMatrixD vecYk(5,1);    // Innovation or measurement residual
  TMatrixD matHkT(5,5);
  TMatrixD matSk(5,5);    // Innovation (or residual) covariance
  TMatrixD matKk(5,5);    // Optimal Kalman gain
  TMatrixD mat1(5,5);     // update covariance matrix
  TMatrixD covXk2(5,5);   // 
  TMatrixD covOut(5,5);
  //
  Double_t *param1=(Double_t*) track1.GetParameter();
  Double_t *covar1=(Double_t*) track1.GetCovariance();
  Double_t *param2=(Double_t*) track2.GetParameter();
  Double_t *covar2=(Double_t*) track2.GetCovariance();
  //
  // copy data to the matrix
  for (Int_t ipar=0; ipar<5; ipar++){
    for (Int_t jpar=0; jpar<5; jpar++){
      covXk(ipar,jpar) = covar1[track1.GetIndex(ipar, jpar)];
      measR(ipar,jpar) = covar2[track2.GetIndex(ipar, jpar)];
      matHk(ipar,jpar)=0;
      mat1(ipar,jpar)=0;
    }
    vecXk(ipar,0) = param1[ipar];
    vecZk(ipar,0) = param2[ipar];
    matHk(ipar,ipar)=1;
    mat1(ipar,ipar)=0;
  }
  //
  //
  //
  //
  //
  vecYk = vecZk-matHk*vecXk;                 // Innovation or measurement residual
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;      // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;              //  Optimal Kalman gain
  vecXk += matKk*vecYk;                      //  updated vector 
  covXk2 = (mat1-(matKk*matHk));
  covOut =  covXk2*covXk; 
  //
  //
  //
  // copy from matrix to parameters
  if (0) {
    vecXk.Print();
    vecZk.Print();
    //
    measR.Print();
    covXk.Print();
    covOut.Print();
    //
    track1.Print();
    track2.Print();
  }

  for (Int_t ipar=0; ipar<5; ipar++){
    param1[ipar]= vecXk(ipar,0) ;
    for (Int_t jpar=0; jpar<5; jpar++){
      covar1[track1.GetIndex(ipar, jpar)]=covOut(ipar,jpar);
    }
  }
}




