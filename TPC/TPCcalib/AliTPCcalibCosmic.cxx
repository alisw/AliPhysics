

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
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVfriendEvent.h"
#include "AliVfriendTrack.h"

#include "AliTracker.h"
#include "AliMagF.h"
#include "AliTPCCalROC.h"
#include "AliTPCParam.h"
#include "AliLog.h"

#include "AliTPCcalibCosmic.h"
#include "TTreeStream.h"
#include "AliTPCTracklet.h"
//#include "AliESDcosmic.h"
#include "AliRecoParam.h"
#include "AliMultiplicity.h"
#include "AliTPCTransform.h"
#include "AliTPCcalibDB.h"
#include "AliTPCseed.h"
#include "AliGRPObject.h"
#include "AliTPCCorrection.h"
#include "AliTPCreco.h"
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
   fCutMinDir(-0.99),   // direction vector products
   fCosmicTree(0)      // tree with cosmic data
{  
  //
  // CONSTRUCTOR - SEE COMMENTS ABOVE
  //
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
   fCutMinDir(-0.99),  // direction vector products
   fCosmicTree(0)      // tree with cosmic data
{  
  //
  // cONSTRUCTOR - SEE COMENTS ABOVE
  //
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
  // destructor
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
  // 9 - is corss indicator
  Int_t ndim=10;
  Double_t xminTrack[10], xmaxTrack[10];
  Int_t binsTrack[10];
  TString axisName[10];
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
  binsTrack[7] =50;
  xminTrack[7] =1; xmaxTrack[7]=1000;  // 
  axisName[7]  ="pt (GeV)";
  //
  binsTrack[8] =18;
  xminTrack[8] =0; xmaxTrack[8]=TMath::Pi();  // 
  axisName[8]  ="alpha";
  //
  binsTrack[9] =3;
  xminTrack[9] =-0.1; xmaxTrack[9]=2.1;  // 
  axisName[9]  ="cross";
  //
  // delta y
  xminTrack[0] =-1; xmaxTrack[0]=1;  // 
  fHistoDelta[0] = new THnSparseS("#Delta_{Y} (cm)","#Delta_{Y} (cm)", ndim, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fHistoPull[0] = new THnSparseS("#Delta_{Y} (unit)","#Delta_{Y} (unit)", ndim, binsTrack,xminTrack, xmaxTrack);
  //
  // delta z
  xminTrack[0] =-1; xmaxTrack[0]=1;  // 
  fHistoDelta[1] = new THnSparseS("#Delta_{Z} (cm)","#Delta_{Z} (cm)", ndim, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fHistoPull[1] = new THnSparseS("#Delta_{Z} (unit)","#Delta_{Z} (unit)", ndim, binsTrack,xminTrack, xmaxTrack);
  //
  // delta P2
  xminTrack[0] =-10; xmaxTrack[0]=10;  // 
  fHistoDelta[2] = new THnSparseS("#Delta_{#phi} (mrad)","#Delta_{#phi} (mrad)", ndim, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fHistoPull[2] = new THnSparseS("#Delta_{#phi} (unit)","#Delta_{#phi} (unit)", ndim, binsTrack,xminTrack, xmaxTrack);
  //
  // delta P3
  xminTrack[0] =-10; xmaxTrack[0]=10;  // 
  fHistoDelta[3] = new THnSparseS("#Delta_{#theta} (mrad)","#Delta_{#theta} (mrad)", ndim, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fHistoPull[3] = new THnSparseS("#Delta_{#theta} (unit)","#Delta_{#theta} (unit)", ndim, binsTrack,xminTrack, xmaxTrack);
  //
  // delta P4
  xminTrack[0] =-0.2; xmaxTrack[0]=0.2;  // 
  fHistoDelta[4] = new THnSparseS("#Delta_{1/pt} (1/GeV)","#Delta_{1/pt} (1/GeV)", ndim, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fHistoPull[4] = new THnSparseS("#Delta_{1/pt} (unit)","#Delta_{1/pt} (unit)", ndim, binsTrack,xminTrack, xmaxTrack);
  
  //
  // delta Pt
  xminTrack[0] =-0.5; xmaxTrack[0]=0.5;  // 
  fHistoDelta[5] = new THnSparseS("#Delta_{pt}/p_{t}","#Delta_{pt}/p_{t}", ndim, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fHistoPull[5] = new THnSparseS("#Delta_{pt}/p_{t} (unit)","#Delta_{pt}/p_{t} (unit)", ndim, binsTrack,xminTrack, xmaxTrack);
  //

  for (Int_t idedx=0;idedx<4;idedx++){
    xminTrack[0] =0.5; xmaxTrack[0]=1.5;  // 
    binsTrack[1] =40;
    xminTrack[1] =10; xmaxTrack[1]=160;

    fHistodEdxMax[idedx] = new THnSparseS(Form("dEdx_{MaxUp}/dEdx_{MaxDown}_Pad%d",idedx),
					  Form("dEdx_{MaxUp}/dEdx_{MaxDown}_Pad%d",idedx), 
					  ndim, binsTrack,xminTrack, xmaxTrack);
    fHistodEdxTot[idedx] = new THnSparseS(Form("dEdx_{TotUp}/dEdx_{TotDown}_Pad%d",idedx),
					  Form("dEdx_{TotUp}/dEdx_{TotDown}_Pad%d",idedx), 
					  ndim, binsTrack,xminTrack, xmaxTrack);
  }
  


  for (Int_t ivar=0;ivar<6;ivar++){
    for (Int_t ivar2=0;ivar2<ndim;ivar2++){      
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
  // merge the content of the cosmic componentnts
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
  if (cosmic->fCosmicTree){
    if (!fCosmicTree) {
      fCosmicTree = new TTree("pairs","pairs");
      fCosmicTree->SetDirectory(0);
    }
    AliTPCcalibCosmic::AddTree(fCosmicTree,cosmic->fCosmicTree);
  }
}




void AliTPCcalibCosmic::Process(AliVEvent *event) {
  //
  // Process of the ESD event  - fill calibration components
  //
  if (!event) {
    Printf("ERROR: event not available");
    return;
  }  
   
  //
  //Int_t isOK=kTRUE;
  // COSMIC not signed properly
  //  UInt_t specie = event->GetEventSpecie();  // select only cosmic events
  //if (specie==AliRecoParam::kCosmic || specie==AliRecoParam::kCalib) {
  //  isOK = kTRUE;
  //}
  //if (!isOK) return;
  // Work around
  FindCosmicPairs(event);
  //const AliMultiplicity *multiplicity = event->GetMultiplicity();
  //  Int_t ntracklets = multiplicity->GetNumberOfTracklets();
  //if (ntracklets>6) return; // filter out "normal" event with high multiplicity
  //const TString &trigger = event->GetFiredTriggerClasses();
  //if (trigger.Contains("C0OB0")==0) return;


  FindPairs(event); // nearly everything takes place in find pairs...

  if (GetDebugLevel()>20) printf("Hallo world: Im here and processing an event\n");
  Int_t ntracks=event->GetNumberOfTracks(); 
  fHistNTracks->Fill(ntracks);
  
}


void AliTPCcalibCosmic::FillHistoPerformance(const AliExternalTrackParam *par0, const AliExternalTrackParam *par1, const AliExternalTrackParam *inner0, const AliExternalTrackParam */*inner1*/, AliTPCseed *seed0,  AliTPCseed *seed1, const AliExternalTrackParam *param0Combined , Int_t cross){
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
  Double_t x[10];
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
  x[9] = cross;
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
    x[0]= delta[ivar];    // Modifiation 10.10 use not normalized deltas
    if (ivar==2 || ivar ==3) x[0]*=1000;  // angles in radian
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

void AliTPCcalibCosmic::FindPairs(const AliVEvent *event){
  //
  // Find cosmic pairs
  // 
  // Track0 is choosen in upper TPC part
  // Track1 is choosen in lower TPC part
  //
  if (GetDebugLevel()>20) printf("Hallo world: Im here\n");
  AliVfriendEvent *friendEvent=event->FindFriend();
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
   AliVTrack *track = event->GetVTrack(i);
   if(!track) continue;
   fClusters->Fill(track->GetTPCNcls()); 
  
   AliExternalTrackParam trckIn;
   if( (track->GetTrackParamIp(trckIn)) < 0 ) continue;
   AliExternalTrackParam * trackIn = &trckIn;

   AliExternalTrackParam trckOut;
   if ( (track->GetTrackParamOp(trckOut)) < 0 ) continue;
   AliExternalTrackParam * trackOut = &trckOut;

   if (ntracks>4 && TMath::Abs(trackIn->GetTgl())<0.0015) continue;  // filter laser 

   const AliVfriendTrack *friendTrack = friendEvent->GetTrack(i);
   if (!friendTrack) continue;
   AliTPCseed *seed = new AliTPCseed();   
   if (friendTrack->GetTPCseed(*seed)==0) {
     tpcSeeds.AddAt(seed,i);
   } 
   else {
     delete seed;
   }

   Double_t meanP = 0.5*(trackIn->GetP() + trackOut->GetP());
   if (seed && track->GetTPCNcls() > 80 + 60/(1+TMath::Exp(-meanP+5))) {
     fDeDx->Fill(meanP, seed->CookdEdxNorm(0.0,0.45,0,0,kMaxRow));
     //
     if (meanP > 0.4 && meanP < 0.45) fDeDxMIP->Fill(seed->CookdEdxNorm(0.0,0.45,0,0,kMaxRow));
     //
     // if (GetDebugLevel()>0&&meanP>0.2&&seed->CookdEdxNorm(0.0,0.45,0,0,kMaxRow)>300) {
//        //TFile *curfile = AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile();
//        //if (curfile) printf(">>> p+ in file: %s \t event: %i \t Number of ESD tracks: %i \n", curfile->GetName(), (int)event->GetEventNumberInFile(), (int)ntracks);
//        // if (track->GetOuterParam()->GetAlpha()<0) cout << " Polartiy: " << track->GetSign() << endl;
//      }

   }

  }

  if (ntracks<2) return;
  //
  // Find pairs
  //
  for (Int_t i=0;i<ntracks;++i) {
      AliVTrack *track0 = event->GetVTrack(i);
    // track0 - choosen upper part
    if (!track0) continue;

    AliExternalTrackParam trk0Out;
    if ( (track0->GetTrackParamOp(trk0Out)) < 0) continue;
    if (trk0Out.GetAlpha()<0) continue;

    Double_t dir0[3];
    track0->GetDirection(dir0);    
    for (Int_t j=0;j<ntracks;++j) {
      if (i==j) continue;
      AliVTrack *track1 = event->GetVTrack(j);
      //track 1 lower part
      if (!track1) continue;

      AliExternalTrackParam trk1Out;
      if ( (track1->GetTrackParamOp(trk1Out)) < 0 ) continue;
      if (trk1Out.GetAlpha()>0) continue;
      //
      Double_t dir1[3];
      track1->GetDirection(dir1);
      
      AliTPCseed * seed0 = (AliTPCseed*) tpcSeeds.At(i);
      AliTPCseed * seed1 = (AliTPCseed*) tpcSeeds.At(j);
      if (! seed0) continue;
      if (! seed1) continue;
      Float_t dedx0 = seed0->CookdEdxNorm(0.05,0.55,0,0,kMaxRow);
      Float_t dedx1 = seed1->CookdEdxNorm(0.05,0.55,0,0,kMaxRow);
      //
      Float_t dedx0I = seed0->CookdEdxNorm(0.05,0.55,0,0,63);
      Float_t dedx1I = seed1->CookdEdxNorm(0.05,0.55,0,0,63);
      //
      Float_t dedx0O = seed0->CookdEdxNorm(0.05,0.55,0,64,kMaxRow);
      Float_t dedx1O = seed1->CookdEdxNorm(0.05,0.55,0,64,kMaxRow);
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
      AliExternalTrackParam param0;
      track0->GetTrackParam(param0);

      AliExternalTrackParam param1;
      track1->GetTrackParam(param1);
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
      if (isPair &&param0.Pt()>1) {
	const TString &trigger = event->GetFiredTriggerClasses(); 
	UInt_t specie = event->GetEventSpecie();
	printf("COSMIC ?\t%s\t%d\t%f\t%f\n", trigger.Data(),specie, param0.GetZ(), param1.GetZ()); 
      }
      //
      // combined track params 
      //
      AliExternalTrackParam *par0U=MakeCombinedTrack(&param0,&param1);
      AliExternalTrackParam *par1U=MakeCombinedTrack(&param1,&param0);
      
      //
      if (fStreamLevel>0){
	TTreeSRedirector * cstream =  GetDebugStreamer();
	//printf("My stream=%p\n",(void*)cstream);
    AliExternalTrackParam trck0ip;
    track0->GetTrackParamIp(trck0ip);
    AliExternalTrackParam *ip0 = &trck0ip;

    AliExternalTrackParam trck1ip;
    track1->GetTrackParamIp(trck1ip);
    AliExternalTrackParam *ip1 = &trck1ip;

    AliExternalTrackParam trck0op;
    track0->GetTrackParamOp(trck0op);
    AliExternalTrackParam *op0 = &trck0op;

    AliExternalTrackParam trck1op;
    track1->GetTrackParamOp(trck1op);
    AliExternalTrackParam *op1 = &trck1op;

	Bool_t isCrossI = ip0->GetZ()*ip1->GetZ()<0;
	Bool_t isCrossO = op0->GetZ()*op1->GetZ()<0;
	Double_t alpha0 = TMath::ATan2(dir0[1],dir0[0]);
	Double_t alpha1 = TMath::ATan2(dir1[1],dir1[0]);
	//
	//
	//
	Int_t cross =0;  // 0 no cross, 2 cross on both sides
	if (isCrossI) cross+=1; 
	if (isCrossO) cross+=1; 
	FillHistoPerformance(&param0, &param1, ip0, ip1, seed0, seed1,par0U, cross);
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
  tpcSeeds.Delete();
}    




void  AliTPCcalibCosmic::FillAcordeHist(AliVTrack *upperTrack) {

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
  //
  // component merging
  //

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


Bool_t  AliTPCcalibCosmic::IsPair(AliExternalTrackParam *tr0, AliExternalTrackParam *tr1) const{
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
  //
  // Calculate the MIP value - gaussian fit used
  //

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
  delete [] newBins;

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
  delete [] newBins;
  
}


AliExternalTrackParam *AliTPCcalibCosmic::MakeTrack(const AliExternalTrackParam *track0, const AliExternalTrackParam *track1){
  //
  // Make a atrack using the kalman update of track0 and track1
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



void AliTPCcalibCosmic::FindCosmicPairs(const AliVEvent *event) {
  //
  // find cosmic pairs trigger by random trigger
  //
  //
  AliESDVertex vtxSPD;
  event->GetPrimaryVertexSPD(vtxSPD);
  AliESDVertex *vertexSPD=&vtxSPD;

  AliESDVertex vtxTPC;
  event->GetPrimaryVertexTPC(vtxTPC);
  AliESDVertex *vertexTPC=&vtxTPC;

  AliVfriendEvent *friendEvent=event->FindFriend();
  const Double_t kMinPt=1;
  const Double_t kMinPtMax=0.8;
  const Double_t kMinNcl=50;
  const Double_t kMaxDelta[5]={2,600,0.02,0.02,0.1};
  Int_t ntracks=event->GetNumberOfTracks(); 
  //  Float_t dcaTPC[2]={0,0};
  // Float_t covTPC[3]={0,0,0};

  UInt_t specie = event->GetEventSpecie();  // skip laser events
  if (specie==AliRecoParam::kCalib) return;
  


  for (Int_t itrack0=0;itrack0<ntracks;itrack0++) {
    AliVTrack *track0 = event->GetVTrack(itrack0);
    if (!track0) continue;
    if (!track0->IsOn(AliVTrack::kTPCrefit)) continue;

    if (TMath::Abs(AliTracker::GetBz())>1&&track0->Pt()<kMinPt) continue;
    if (track0->GetTPCncls()<kMinNcl) continue;
    if (TMath::Abs(track0->GetY())<kMaxDelta[0]) continue; 
    if (track0->GetKinkIndex(0)>0) continue;

    AliExternalTrackParam trkprm0;
    track0->GetTrackParam(trkprm0);
    const Double_t * par0=trkprm0.GetParameter(); //track param at the DCA

    //rm primaries
    //
    //track0->GetImpactParametersTPC(dcaTPC,covTPC);
    //if (TMath::Abs(dcaTPC[0])<kMaxDelta[0]) continue;
    //if (TMath::Abs(dcaTPC[1])<kMaxDelta[0]*2) continue;
    //    const AliExternalTrackParam * trackIn0 = track0->GetInnerParam();
    for (Int_t itrack1=itrack0+1;itrack1<ntracks;itrack1++) {
      AliVTrack *track1 = event->GetVTrack(itrack1);
      if (!track1) continue;  
      if (!track1->IsOn(AliVTrack::kTPCrefit)) continue;
      if (track1->GetKinkIndex(0)>0) continue;
      if (TMath::Abs(AliTracker::GetBz())>1&&track1->Pt()<kMinPt) continue;
      if (track1->GetTPCncls()<kMinNcl) continue;
      if (TMath::Abs(AliTracker::GetBz())>1&&TMath::Max(track1->Pt(), track0->Pt())<kMinPtMax) continue;
      if (TMath::Abs(track1->GetY())<kMaxDelta[0]) continue;
      //track1->GetImpactParametersTPC(dcaTPC,covTPC);
      //      if (TMath::Abs(dcaTPC[0])<kMaxDelta[0]) continue;
      //if (TMath::Abs(dcaTPC[1])<kMaxDelta[0]*2) continue;
      //
      AliExternalTrackParam trkprm1;
      track1->GetTrackParam(trkprm1);
      const Double_t* par1=trkprm1.GetParameter(); //track param at the DCA
      //
      Bool_t isPair=kTRUE;
      for (Int_t ipar=0; ipar<5; ipar++){
	if (ipar==4&&TMath::Abs(AliTracker::GetBz())<1) continue; // 1/pt not defined for B field off
	if (TMath::Abs(TMath::Abs(par0[ipar])-TMath::Abs(par1[ipar]))>kMaxDelta[ipar]) isPair=kFALSE;
      }
      if (!isPair) continue;
      if (TMath::Abs(TMath::Abs(track0->GetAlpha()-track1->GetAlpha())-TMath::Pi())>kMaxDelta[2]) isPair=kFALSE;
      //delta with correct sign
      /*
	TCut cut0="abs(t1.fP[0]+t0.fP[0])<2"
	TCut cut3="abs(t1.fP[3]+t0.fP[3])<0.02"
	TCut cut4="abs(t1.fP[4]+t0.fP[4])<0.2"
      */
      if  (TMath::Abs(par0[0]+par1[0])>kMaxDelta[0]) isPair=kFALSE; //delta y   opposite sign
      if  (TMath::Abs(par0[3]+par1[3])>kMaxDelta[3]) isPair=kFALSE; //delta tgl opposite sign
      if  (TMath::Abs(AliTracker::GetBz())>1 && TMath::Abs(par0[4]+par1[4])>kMaxDelta[4]) isPair=kFALSE; //delta 1/pt opposite sign
      if (!isPair) continue;
      TString filename(AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()->GetName());
      Int_t eventNumber = event->GetEventNumberInFile(); 
      Bool_t hasFriend=(friendEvent) ? (friendEvent->GetTrack(itrack0)!=0):0;
      Bool_t hasITS=(track0->GetNcls(0)+track1->GetNcls(0)>4);
      printf("DUMPHPTCosmic:%s|%f|%d|%d|%d\n",filename.Data(),(TMath::Min(track0->Pt(),track1->Pt())), eventNumber,hasFriend,hasITS);


      //      const AliExternalTrackParam * trackIn1 = track1->GetInnerParam();      
      //
      //       
      TTreeSRedirector * pcstream =  GetDebugStreamer();
      Int_t ntracksSPD = vertexSPD->GetNContributors();
      Int_t ntracksTPC = vertexTPC->GetNContributors();
      //
      AliVfriendTrack *friendTrack0 = const_cast<AliVfriendTrack*>(friendEvent->GetTrack(itrack0));
      if (!friendTrack0) continue;
      AliVfriendTrack *friendTrack1 = const_cast<AliVfriendTrack*>(friendEvent->GetTrack(itrack1));
      if (!friendTrack1) continue;
      AliTPCseed *seed0 = 0;   
      AliTPCseed *seed1 = 0;
      AliTPCseed tpcSeed0;
      AliTPCseed tpcSeed1;
      //
      if (friendTrack0->GetTPCseed(tpcSeed0)==0) seed0=&tpcSeed0;
      if (friendTrack1->GetTPCseed(tpcSeed1)==0) seed1=&tpcSeed1;

      //
      if (pcstream){
	(*pcstream)<<"pairs"<<
	  "run="<<fRun<<              //  run number
	  "event="<<fEvent<<          //  event number
	  "time="<<fTime<<            //  time stamp of event
	  "trigger="<<fTrigger<<      //  trigger
	  "triggerClass="<<&fTriggerClass<<      //  trigger
	  "mag="<<fMagF<<             //  magnetic field
	  //
	  "nSPD="<<ntracksSPD<<
	  "nTPC="<<ntracksTPC<<
	  "vSPD.="<<vertexSPD<<         //primary vertex -SPD
	  "vTPC.="<<vertexTPC<<         //primary vertex -TPC
	  "t0.="<<track0<<              //track0
	  "t1.="<<track1<<              //track1
	  "ft0.="<<friendTrack0<<       //track0
	  "ft1.="<<friendTrack1<<       //track1
 	  "s0.="<<seed0<<               //track0
 	  "s1.="<<seed1<<               //track1
	  "\n";      
      }

      if (!fCosmicTree) {
	fCosmicTree = new TTree("pairs","pairs");
	fCosmicTree->SetDirectory(0);
      }
      if (fCosmicTree->GetEntries()==0){
	//
	fCosmicTree->SetDirectory(0);
    fCosmicTree->Branch("t0.",&track0);
    fCosmicTree->Branch("t1.",&track1);
    fCosmicTree->Branch("ft0.",&friendTrack0);
    fCosmicTree->Branch("ft1.",&friendTrack1);
      }else{
    fCosmicTree->SetBranchAddress("t0.",&track0);
    fCosmicTree->SetBranchAddress("t1.",&track1);
    fCosmicTree->SetBranchAddress("ft0.",&friendTrack0);
    fCosmicTree->SetBranchAddress("ft1.",&friendTrack1);
      }
      fCosmicTree->Fill();
    }
  }
}


void  AliTPCcalibCosmic::Terminate(){
  //
  // copy the cosmic tree to memory resident tree
  //
  static Int_t counter=0;
  printf("AliTPCcalibCosmic::Terminate\t%d\n",counter);
  counter++;
  AliTPCcalibBase::Terminate();
}


void AliTPCcalibCosmic::AddTree(TTree* treeOutput, TTree * treeInput){
  //
  // Add the content of tree: 
  // Notice automatic copy of tree in ROOT does not work for such complicated tree
  //  
  return;
  //if (TMath::Abs(fMagF)<0.1) return; // work around - otherwise crashes 
  AliVTrack *track0 = 0;
  AliVTrack *track1 = 0;
  AliVfriendTrack *ftrack0 = 0;
  AliVfriendTrack *ftrack1 = 0;
  treeInput->SetBranchAddress("t0.",&track0);	
  treeInput->SetBranchAddress("t1.",&track1);
  treeInput->SetBranchAddress("ft0.",&ftrack0);	
  treeInput->SetBranchAddress("ft1.",&ftrack1);
  treeOutput->SetDirectory(0);
  //
  Int_t entries= treeInput->GetEntries();
  Int_t step=1+Int_t(TMath::Log(1+treeOutput->GetEntries()/10.));
  for (Int_t i=0; i<entries; i+=step){
    treeInput->SetBranchAddress("t0.",&track0);	
    treeInput->SetBranchAddress("t1.",&track1);
    treeInput->SetBranchAddress("ft0.",&ftrack0);	
    treeInput->SetBranchAddress("ft1.",&ftrack1);
    treeInput->GetEntry(i);
    if (!track0) continue;
    if (!track1) continue;
    if (!ftrack0) continue;
    if (!ftrack1) continue;
    if (track0->GetTPCncls()<=0) continue;
    if (track1->GetTPCncls()<=0) continue;

    AliExternalTrackParam trck0Ip;
    if ( (track0->GetTrackParamIp(trck0Ip)) < 0 ) continue;
    AliExternalTrackParam trck1Ip;
    if ( (track1->GetTrackParamIp(trck1Ip)) < 0 ) continue;
    AliExternalTrackParam trck0TPCIn;
    if ( (track0->GetTrackParamTPCInner(trck0TPCIn)) < 0 ) continue;
    AliExternalTrackParam trck1TPCIn;
    if ( (track1->GetTrackParamTPCInner(trck1TPCIn)) < 0 ) continue;
    //track0
    treeOutput->SetBranchAddress("t0.",&track0);	
    treeOutput->SetBranchAddress("t1.",&track1);
    treeOutput->SetBranchAddress("ft0.",&ftrack0);	
    treeOutput->SetBranchAddress("ft1.",&ftrack1);    
    treeOutput->Fill();
    delete track0;
    delete track1;
    delete ftrack0;
    delete ftrack1;
    track0=0;
    track1=0;
    ftrack0=0;
    ftrack1=0;
  }
}



void AliTPCcalibCosmic::MakeFitTree(TTree * treeInput, TTreeSRedirector *pcstream, const TObjArray * corrArray, Int_t step, Int_t run){
  //
  // Make fit tree
  // refit the tracks with original points + corrected points for each correction
  // Input:
  //   treeInput - tree with cosmic tracks
  //   pcstream  - debug output

  // Algorithm:
  // Loop over pair of cosmic tracks:
  //   1. Find trigger offset between cosmic event and "physic" trigger
  //   2. Refit tracks with current transformation
  //   3. Refit tracks using additional "primitive" distortion on top
  // Best correction estimated as linear combination of corrections 
  // minimizing the observed distortions
  // Observed distortions - matching in the y,z, snp, theta and 1/pt
  //
  const Double_t kResetCov=20.;
  const Double_t kMaxDelta[5]={1,40,0.03,0.01,0.2};
  const Double_t kSigma=2.;    
  const Double_t kMaxTime=1050;
  const Double_t kMaxSnp=0.8;
  Int_t ncorr=corrArray->GetEntries();
  AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform() ;
  AliTPCParam     *param     = AliTPCcalibDB::Instance()->GetParameters();
  AliGRPObject*  grp = AliTPCcalibDB::Instance()->GetGRP(run);
  Double_t time=0.5*(grp->GetTimeStart() +grp->GetTimeEnd()); 
  transform->SetCurrentRun(run);
  transform->SetCurrentTimeStamp(TMath::Nint(time));
  Double_t covar[15];
  for (Int_t i=0;i<15;i++) covar[i]=0;
  covar[0]=kSigma*kSigma;
  covar[2]=kSigma*kSigma;
  covar[5]=kSigma*kSigma/Float_t(150*150);
  covar[9]=kSigma*kSigma/Float_t(150*150);
  covar[14]=0.2*0.2;
  Double_t *distortions = new Double_t[ncorr+1];

  AliVTrack *track0 = 0;
  AliVTrack *track1 = 0;
  AliVfriendTrack *ftrack0 = 0;
  AliVfriendTrack *ftrack1 = 0;
  treeInput->SetBranchAddress("t0.",&track0);	
  treeInput->SetBranchAddress("t1.",&track1);
  treeInput->SetBranchAddress("ft0.",&ftrack0);	
  treeInput->SetBranchAddress("ft1.",&ftrack1);
  Int_t entries= treeInput->GetEntries();
  for (Int_t i=0; i<entries; i+=step){    
    treeInput->GetEntry(i);
    if (i%20==0) printf("%d\n",i);

    AliExternalTrackParam ftrck0TPCOut;
    if ( (ftrack0->GetTrackParamTPCOut(ftrck0TPCOut)) < 0 ) continue;
    AliExternalTrackParam ftrck1TPCOut;
    if ( (ftrack1->GetTrackParamTPCOut(ftrck1TPCOut)) < 0 ) continue;
    AliTPCseed *seed0=0;
    AliTPCseed *seed1=0;
    AliTPCseed tpcSeed0;
    AliTPCseed tpcSeed1;
    if (ftrack0->GetTPCseed(tpcSeed0)==0) seed0=&tpcSeed0;
    if (ftrack1->GetTPCseed(tpcSeed1)==1) seed1=&tpcSeed1;
    if (!seed0) continue;
    if (!seed1) continue;
    if (TMath::Abs(seed0->GetSnp())>kMaxSnp) continue;
    if (TMath::Abs(seed1->GetSnp())>kMaxSnp) continue;
    //
    //
    Int_t nclA0=0, nclC0=0;     // number of clusters
    Int_t nclA1=0, nclC1=0;     // number of clusters
    Int_t ncl0=0,ncl1=0;
    Double_t rmin0=300, rmax0=-300;  // variables to estimate the time0 - trigger offset
    Double_t rmin1=300, rmax1=-300;
    Double_t tmin0=2000, tmax0=-2000;
    Double_t tmin1=2000, tmax1=-2000;
    //
    //
    // calculate trigger offset usig "missing clusters"
    for (Int_t irow=0; irow<kMaxRow; irow++){
      AliTPCclusterMI *cluster0=seed0->GetClusterPointer(irow);
      if (cluster0 &&cluster0->GetX()>10){
	if (cluster0->GetX()<rmin0) { rmin0=cluster0->GetX(); tmin0=cluster0->GetTimeBin();}
	if (cluster0->GetX()>rmax0) { rmax0=cluster0->GetX(); tmax0=cluster0->GetTimeBin();}
	ncl0++;
	if (cluster0->GetDetector()%36<18) nclA0++;
	if (cluster0->GetDetector()%36>=18) nclC0++;
      }  
      AliTPCclusterMI *cluster1=seed1->GetClusterPointer(irow);
      if (cluster1&&cluster1->GetX()>10){
	if (cluster1->GetX()<rmin1) { rmin1=cluster1->GetX();  tmin1=cluster1->GetTimeBin();}
	if (cluster1->GetX()>rmax1) { rmax1=cluster1->GetX(); tmax1=cluster1->GetTimeBin();}
	ncl1++;
	if (cluster1->GetDetector()%36<18) nclA1++;
	if (cluster1->GetDetector()%36>=18) nclC1++;
      }
    }
    Int_t cosmicType=0;  // types of cosmic topology
    if ((nclA0>nclC0) && (nclA1>nclC1)) cosmicType=0; // AA side
    if ((nclA0<nclC0) && (nclA1<nclC1)) cosmicType=1; // CC side
    if ((nclA0>nclC0) && (nclA1<nclC1)) cosmicType=2; // AC side
    if ((nclA0<nclC0) && (nclA1>nclC1)) cosmicType=3; // CA side
    //if ((nclA0>nclC0) && (nclA1<nclC1)) cosmicType=6; // AC side out of time
    //if ((nclA0>nclC0) && (nclA1<nclC1)) cosmicType=7; // CA side out of time
    //
    Double_t deltaTime = 0;   // correction for trigger not in time - equalizing the time arival
    if ((cosmicType>1)&&TMath::Abs(track1->GetZ()-track0->GetZ())>4){
      cosmicType+=4;
      deltaTime=0.5*(track1->GetZ()-track0->GetZ())/param->GetZWidth();
      if (nclA0>nclC0) deltaTime*=-1; // if A side track
    }
    //
    TVectorD vectorDT(8);
    Int_t crossCounter=0;
    Double_t deltaTimeCross = AliTPCcalibCosmic::GetDeltaTime(rmin0, rmax0, rmin1, rmax1, tmin0, tmax0, tmin1, tmax1, TMath::Abs(track0->GetY()),vectorDT);
    Bool_t isOKTrigger=kTRUE;
    for (Int_t ic=0; ic<6;ic++) {
      if (TMath::Abs(vectorDT[ic])>0) {
	if (vectorDT[ic]+vectorDT[6]<0) isOKTrigger=kFALSE;
	if (vectorDT[ic]+vectorDT[7]>kMaxTime) isOKTrigger=kFALSE;
	if (isOKTrigger){
	  crossCounter++; 
	}
      }
    }
    Double_t deltaTimeCluster=deltaTime;
    if ((cosmicType==0 || cosmicType==1) && crossCounter>0){
      deltaTimeCluster=deltaTimeCross;
      cosmicType+=8;
    }
    if (nclA0*nclC0>0 || nclA1*nclC1>0) cosmicType+=16;  // mixed A side C side - bad for visualization
    //
    // Apply current transformation
    //
    //
    for (Int_t irow=0; irow<kMaxRow; irow++){
      AliTPCclusterMI *cluster0=seed0->GetClusterPointer(irow);
      if (cluster0 &&cluster0->GetX()>10){
	Double_t x0[3]={ static_cast<Double_t>(cluster0->GetRow()),cluster0->GetPad(),cluster0->GetTimeBin()+deltaTimeCluster};
	Int_t index0[1]={cluster0->GetDetector()};
	transform->Transform(x0,index0,0,1);  
	cluster0->SetX(x0[0]);
	cluster0->SetY(x0[1]);
	cluster0->SetZ(x0[2]);
	//
      }
      AliTPCclusterMI *cluster1=seed1->GetClusterPointer(irow);
      if (cluster1&&cluster1->GetX()>10){
	Double_t x1[3]={ static_cast<Double_t>(cluster1->GetRow()),cluster1->GetPad(),cluster1->GetTimeBin()+deltaTimeCluster};
	Int_t index1[1]={cluster1->GetDetector()};
	transform->Transform(x1,index1,0,1);  
	cluster1->SetX(x1[0]);
	cluster1->SetY(x1[1]);
	cluster1->SetZ(x1[2]);
      }
    }
    //
    //
    Double_t alpha=track0->GetAlpha();   // rotation frame
    Double_t cos = TMath::Cos(alpha);
    Double_t sin = TMath::Sin(alpha);
    Double_t mass =  TDatabasePDG::Instance()->GetParticle("mu+")->Mass();
    AliExternalTrackParam  btrack0;
    ftrack0->GetTrackParamTPCOut(btrack0);
    AliExternalTrackParam  btrack1;
    ftrack1->GetTrackParamTPCOut(btrack1);
    btrack0.Rotate(alpha);
    btrack1.Rotate(alpha);
    // change the sign for track 1
    Double_t* par1=(Double_t*)btrack0.GetParameter();
    par1[3]*=-1;
    par1[4]*=-1;
    btrack0.AddCovariance(covar);
    btrack1.AddCovariance(covar);
    btrack0.ResetCovariance(kResetCov);
    btrack1.ResetCovariance(kResetCov);
    Bool_t isOK=kTRUE;
    Bool_t isOKT=kTRUE;
    TObjArray tracks0(ncorr+1);
    TObjArray tracks1(ncorr+1);
    //    
    Double_t dEdx0Tot=seed0->CookdEdxAnalytical(0.02,0.6,kTRUE);
    Double_t dEdx1Tot=seed1->CookdEdxAnalytical(0.02,0.6,kTRUE);
    Double_t dEdx0Max=seed0->CookdEdxAnalytical(0.02,0.6,kFALSE);
    Double_t dEdx1Max=seed1->CookdEdxAnalytical(0.02,0.6,kFALSE);
    //if (TMath::Abs((dEdx0Max+1)/(dEdx0Tot+1)-1.)>0.1) isOK=kFALSE;
    //if (TMath::Abs((dEdx1Max+1)/(dEdx1Tot+1)-1.)>0.1) isOK=kFALSE;
    ncl0=0; ncl1=0;
    for (Int_t icorr=-1; icorr<ncorr; icorr++){
      AliExternalTrackParam  rtrack0=btrack0;
      AliExternalTrackParam  rtrack1=btrack1;
      AliTPCCorrection *corr = 0;
      if (icorr>=0) corr = (AliTPCCorrection*)corrArray->At(icorr);
      //
      for (Int_t irow=kMaxRow; irow--;){ 
	AliTPCclusterMI *cluster=seed0->GetClusterPointer(irow);
	if (!cluster) continue;
	if (!isOKT) break;
	Double_t rD[3]={cluster->GetX(),cluster->GetY(),cluster->GetZ()};
	transform->RotatedGlobal2Global(cluster->GetDetector()%36,rD);  // transform to global
	Float_t  r[3]={static_cast<Float_t>(rD[0]),static_cast<Float_t>(rD[1]),static_cast<Float_t>(rD[2])};
	if (corr){
	  corr->DistortPoint(r, cluster->GetDetector());
	}
	//
	Double_t cov[3]={0.01,0.,0.01}; 
	Double_t lx =cos*r[0]+sin*r[1];      
	Double_t ly =-sin*r[0]+cos*r[1];
	rD[1]=ly; rD[0]=lx; rD[2]=r[2];  //transform to track local
	if (!AliTracker::PropagateTrackToBxByBz(&rtrack0, lx,mass,1.,kFALSE)) isOKT=kFALSE;;
	if (!rtrack0.Update(&rD[1],cov)) isOKT =kFALSE;
	if (icorr<0) ncl0++;
      }
      //
      for (Int_t irow=kMaxRow; irow--;){ 
	AliTPCclusterMI *cluster=seed1->GetClusterPointer(irow);
	if (!cluster) continue;
	if (!isOKT) break;
	Double_t rD[3]={cluster->GetX(),cluster->GetY(),cluster->GetZ()};
	transform->RotatedGlobal2Global(cluster->GetDetector()%36,rD);
	Float_t  r[3]={static_cast<Float_t>(rD[0]),static_cast<Float_t>(rD[1]),static_cast<Float_t>(rD[2])};
	if (corr){
	  corr->DistortPoint(r, cluster->GetDetector());
	}
	//
	Double_t cov[3]={0.01,0.,0.01}; 
	Double_t lx =cos*r[0]+sin*r[1];      
	Double_t ly =-sin*r[0]+cos*r[1];
	rD[1]=ly; rD[0]=lx; rD[2]=r[2];
	if (!AliTracker::PropagateTrackToBxByBz(&rtrack1, lx,mass,1.,kFALSE)) isOKT=kFALSE;
	if (!rtrack1.Update(&rD[1],cov)) isOKT=kFALSE;
	if (icorr<0) ncl1++;
      }
      if (!AliTracker::PropagateTrackToBxByBz(&rtrack0, 0,mass,10.,kFALSE)) isOKT=kFALSE;
      if (!AliTracker::PropagateTrackToBxByBz(&rtrack1, 0,mass,10.,kFALSE)) isOKT=kFALSE;
      if (!AliTracker::PropagateTrackToBxByBz(&rtrack0, 0,mass,1.,kFALSE))  isOKT=kFALSE;
      if (!AliTracker::PropagateTrackToBxByBz(&rtrack1, 0,mass,1.,kFALSE))  isOKT=kFALSE;
      const Double_t *param0=rtrack0.GetParameter();
      const Double_t *param1=rtrack1.GetParameter();
      for (Int_t ipar=0; ipar<4; ipar++){
	if (TMath::Abs(param1[ipar]-param0[ipar])>kMaxDelta[ipar]) isOK=kFALSE;
      }
      tracks0.AddAt(rtrack0.Clone(), icorr+1);
      tracks1.AddAt(rtrack1.Clone(), icorr+1);
      AliExternalTrackParam out0;
      ftrack0->GetTrackParamTPCOut(out0);
      AliExternalTrackParam out1;
      ftrack1->GetTrackParamTPCOut(out1);
      Int_t nentries=TMath::Min(ncl0,ncl1);

      if (icorr<0) {
	(*pcstream)<<"cosmic"<<
	  "isOK="<<isOK<<              // correct all propagation update and also residuals
	  "isOKT="<<isOKT<<            // correct all propagation update
	  "isOKTrigger="<<isOKTrigger<< // correct? estimate of trigger offset
	  "id="<<cosmicType<<
	  //
	  //
	  "cross="<<crossCounter<<
	  "vDT.="<<&vectorDT<<
	  //
	  "dTime="<<deltaTime<<        // delta time using the A-c side cross
	  "dTimeCross="<<deltaTimeCross<< // delta time using missing clusters
	  //
	  "dEdx0Max="<<dEdx0Max<<
	  "dEdx0Tot="<<dEdx0Tot<<
	  "dEdx1Max="<<dEdx1Max<<
	  "dEdx1Tot="<<dEdx1Tot<<
	  //
	  "track0.="<<track0<<         // original track 0
	  "track1.="<<track1<<         // original track 1
	  "out0.="<<&out0<<             // outer track  0
	  "out1.="<<&out1<<             // outer track  1
	  "rtrack0.="<<&rtrack0<<      // refitted track with current transform
	  "rtrack1.="<<&rtrack1<<     //	  
	  "ncl0="<<ncl0<<
	  "ncl1="<<ncl1<<
	  "entries="<<nentries<<       // number of clusters
	  "\n";
      }
    }
    //

    if (isOK){        
      Int_t nentries=TMath::Min(ncl0,ncl1);    
      for (Int_t ipar=0; ipar<5; ipar++){
	for (Int_t icorr=-1; icorr<ncorr; icorr++){
	  AliTPCCorrection *corr = 0;
	  if (icorr>=0) corr = (AliTPCCorrection*)corrArray->At(icorr);
	  //
	  AliExternalTrackParam *param0=(AliExternalTrackParam *) tracks0.At(icorr+1);
	  AliExternalTrackParam *param1=(AliExternalTrackParam *) tracks1.At(icorr+1);
	  distortions[icorr+1]=param1->GetParameter()[ipar]-param0->GetParameter()[ipar];
	  if (icorr>=0){
	    distortions[icorr+1]-=distortions[0];
	  }
	  //
	  if (icorr<0){
	    Double_t bz=AliTrackerBase::GetBz();
	    Double_t gxyz[3];
	    param0->GetXYZ(gxyz);
	    Int_t dtype=20;
	    Double_t theta=param0->GetParameter()[3];
	    Double_t phi = alpha;
        AliExternalTrackParam trk0Ip;
        track0->GetTrackParamIp(trk0Ip);
        Double_t snp = trk0Ip.GetSnp();
	    Double_t mean= distortions[0];
	    Int_t index = param0->GetIndex(ipar,ipar);
	    Double_t rms=TMath::Sqrt(param1->GetCovariance()[index]+param1->GetCovariance()[index]);
	    if (crossCounter<1) rms*=1;
	    Double_t sector=9*phi/TMath::Pi();
	    Double_t dsec   = sector-TMath::Nint(sector+0.5);
	    Double_t gx=gxyz[0],gy=gxyz[1],gz=gxyz[2];
	    Double_t refX=TMath::Sqrt(gx*gx+gy*gy);
	    Double_t dRrec=0;
	    //	    Double_t pt=(param0->GetSigned1Pt()+param1->GetSigned1Pt())*0.5;
	    Double_t pt=(param0->GetSigned1Pt()+param1->GetSigned1Pt())*0.5;

	    (*pcstream)<<"fit"<<  // dump valus for fit
	      "run="<<run<<       //run number
	      "bz="<<bz<<         // magnetic filed used
	      "dtype="<<dtype<<   // detector match type 20
	      "ptype="<<ipar<<   // parameter type
	      "theta="<<theta<<   // theta
	      "phi="<<phi<<       // phi 
	      "snp="<<snp<<       // snp
	      "mean="<<mean<<     // mean dist value
	      "rms="<<rms<<       // rms
	      "sector="<<sector<<
	      "dsec="<<dsec<<
	      //
	      "refX="<<refX<<      // reference radius
	      "gx="<<gx<<         // global position
	      "gy="<<gy<<         // global position
	      "gz="<<gz<<         // global position
	      "dRrec="<<dRrec<<      // delta Radius in reconstruction
	      "pt="<<pt<<            //1/pt
	      "id="<<cosmicType<<     //type of the cosmic used      
	      "entries="<<nentries;// number of entries in bin
	    (*pcstream)<<"cosmicDebug"<<
	      "p0.="<<param0<<    // dump distorted track 0
	      "p1.="<<param1;    // dump distorted track 1
	  }
	  if (icorr>=0){
	    (*pcstream)<<"fit"<<
	      Form("%s=",corr->GetName())<<distortions[icorr+1];    // dump correction value
	    (*pcstream)<<"cosmicDebug"<<
	      Form("%s=",corr->GetName())<<distortions[icorr+1]<<    // dump correction value
	      Form("%sp0.=",corr->GetName())<<param0<<    // dump distorted track 0
	      Form("%sp1.=",corr->GetName())<<param1;    // dump distorted track 1
	  }
	} //loop corrections      
	(*pcstream)<<"fit"<<"isOK="<<isOK<<"\n";
	(*pcstream)<<"cosmicDebug"<<"isOK="<<isOK<<"\n";
      } //loop over parameters
    } // dump results
  }//loop tracks
  delete [] distortions;
}



Double_t AliTPCcalibCosmic::GetDeltaTime(Double_t rmin0, Double_t rmax0, Double_t rmin1, Double_t rmax1, Double_t tmin0, Double_t tmax0, Double_t tmin1, Double_t tmax1, Double_t dcaR, TVectorD &vectorDT)
{
  //
  // Estimate trigger offset between random cosmic event and "physics" trigger
  // Efficiency about 50 % of cases:
  // Cases:
  // 0. Tracks crossing A side and C side - no match in z - 30 % of cases
  // 1. Track only on one side and  dissapear at small or at high radiuses - 50 % of cases
  //    1.a) rmax<Rc    && tmax<Tcmax && tmax>tmin    => deltaT=Tcmax-tmax 
  //    1.b) rmin>Rcmin && tmin<Tcmax && tmin>tmax    => deltaT=Tcmax-tmin  
  //                      // case the z matching gives proper time
  //    1.c) rmax<Rc    && tmax>Tcmin && tmax<tmin    => deltaT=-tmax
  //
  // check algorithm:
  //    TCut cutStop = "(min(rmax0,rmax1)<235||abs(rmin0-rmin1)>10)"; // tracks not registered for full time

  // Combinations:
  // 0-1 - forbidden
  // 0-2 - forbidden
  // 0-3 - occur - wrong correlation
  // 1-2 - occur - wrong correlation
  // 1-3 - forbidden
  // 2-3 - occur - small number of outlyers -20%
  // Frequency:
  //   0 - 106
  //   1 - 265
  //   2 - 206
  //   3 - 367
  //
  const Double_t kMaxRCut=235;  // max radius
  const Double_t kMinRCut=TMath::Max(dcaR,90.);  // min radius
  const Double_t kMaxDCut=30;   // max distance for minimal radius
  const Double_t kMinTime=110;
  const Double_t kMaxTime=950;  
  Double_t deltaT=0;
  Int_t counter=0;
  vectorDT[6]=TMath::Min(TMath::Min(tmin0,tmin1),TMath::Min(tmax0,tmax1));
  vectorDT[7]=TMath::Max(TMath::Max(tmin0,tmin1),TMath::Max(tmax0,tmax1));
  if (TMath::Min(rmax0,rmax1)<kMaxRCut){
    // max cross - deltaT>0
    if (rmax0<kMaxRCut && tmax0 <kMaxTime && tmax0>tmin0) vectorDT[0]=kMaxTime-tmax0; // disapear at CE
    if (rmax1<kMaxRCut && tmax1 <kMaxTime && tmax1>tmin1) vectorDT[1]=kMaxTime-tmax1; // disapear at CE
    // min cross - deltaT<0 - OK they are correlated
    if (rmax0<kMaxRCut && tmax0 >kMinTime && tmax0<tmin0) vectorDT[2]=-tmax0;         // disapear at ROC
    if (rmax1<kMaxRCut && tmax1 >kMinTime && tmax1<tmin1) vectorDT[3]=-tmax1;         // disapear at ROC
  }  
  if (rmin0> kMinRCut+kMaxDCut && tmin0 <kMaxTime && tmin0>tmax0) vectorDT[4]=kMaxTime-tmin0;
  if (rmin1> kMinRCut+kMaxDCut && tmin1 <kMaxTime && tmin1>tmax1) vectorDT[5]=kMaxTime-tmin1;
  Bool_t isOK=kTRUE;
  for (Int_t i=0; i<6;i++) {
    if (TMath::Abs(vectorDT[i])>0) {
      counter++; 
      if (vectorDT[i]+vectorDT[6]<0) isOK=kFALSE;
      if (vectorDT[i]+vectorDT[7]>kMaxTime) isOK=kFALSE;
      if (isOK) deltaT=vectorDT[i];
    }
  }
  return deltaT;  
}

