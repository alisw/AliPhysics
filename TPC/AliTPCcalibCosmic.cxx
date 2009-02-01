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
#include "AliExternalComparison.h"
#include "TTreeStream.h"
#include "AliTPCTracklet.h"

ClassImp(AliTPCcalibCosmic)


AliTPCcalibCosmic::AliTPCcalibCosmic() 
  :AliTPCcalibBase(),
   fGainMap(0),
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
}


AliTPCcalibCosmic::AliTPCcalibCosmic(const Text_t *name, const Text_t *title) 
  :AliTPCcalibBase(),
   fGainMap(0),
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

  AliInfo("Non Default Constructor");  
  //
}

AliTPCcalibCosmic::~AliTPCcalibCosmic(){
  //
  //
  //
}





void AliTPCcalibCosmic::Process(AliESDEvent *event) {
  //
  //
  //
  if (!event) {
    Printf("ERROR: ESD not available");
    return;
  }  
  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!ESDfriend) {
   Printf("ERROR: ESDfriend not available");
   return;
  }
   

  FindPairs(event); // nearly everything takes place in find pairs...

  if (GetDebugLevel()>20) printf("Hallo world: Im here and processing an event\n");
  Int_t ntracks=event->GetNumberOfTracks(); 
  fHistNTracks->Fill(ntracks);
  if (ntracks==0) return;

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
  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
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


   AliESDfriendTrack *friendTrack = ESDfriend->GetTrack(i);
   TObject *calibObject;
   AliTPCseed *seed = 0;
   for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
     if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
   }
   if (seed) tpcSeeds.AddAt(seed,i);

   Double_t meanP = 0.5*(trackIn->GetP() + trackOut->GetP());
   if (seed && track->GetTPCNcls() > 80 + 60/(1+TMath::Exp(-meanP+5))) {
     fDeDx->Fill(meanP, seed->CookdEdxNorm(0.0,0.45,0,0,159,fGainMap));
     //
     if (meanP > 0.4 && meanP < 0.45) fDeDxMIP->Fill(seed->CookdEdxNorm(0.0,0.45,0,0,159,fGainMap));
     //
     if (GetDebugLevel()>0&&meanP>0.2&&seed->CookdEdxNorm(0.0,0.45,0,0,159,fGainMap)>300) {
       TFile *curfile = AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile();
       if (curfile) printf(">>> p+ in file: %s \t event: %i \t Number of ESD tracks: %i \n", curfile->GetName(), (int)event->GetEventNumberInFile(), (int)ntracks);
       if (track->GetOuterParam()->GetAlpha()<0) cout << " Polartiy: " << track->GetSign() << endl;
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
      Float_t dedx0 = seed0->CookdEdxNorm(0.05,0.55,0,0,159,fGainMap);
      Float_t dedx1 = seed1->CookdEdxNorm(0.05,0.55,0,0,159,fGainMap);
      //
      Float_t dedx0I = seed0->CookdEdxNorm(0.05,0.55,0,0,63,fGainMap);
      Float_t dedx1I = seed1->CookdEdxNorm(0.05,0.55,0,0,63,fGainMap);
      //
      Float_t dedx0O = seed0->CookdEdxNorm(0.05,0.55,0,64,159,fGainMap);
      Float_t dedx1O = seed1->CookdEdxNorm(0.05,0.55,0,64,159,fGainMap);
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
      AliTracker::PropagateTrackTo(&param0,dmax+1,0.0005,3,kTRUE);
      AliTracker::PropagateTrackTo(&param1,dmax+1,0.0005,3,kTRUE);
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
	if (cstream) {
	  (*cstream) << "Track0" <<
	    "run="<<fRun<<              //  run number
	    "event="<<fEvent<<          //  event number
	    "time="<<fTime<<            //  time stamp of event
	    "trigger="<<fTrigger<<      //  trigger
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
    }
  }  
}    




void  AliTPCcalibCosmic::FillAcordeHist(AliESDtrack *upperTrack) {

  // Pt cut to select straight tracks which can be easily propagated to ACORDE which is outside the magnetic field
  if (upperTrack->Pt() < 10 || upperTrack->GetTPCNcls() < 80) return;
    
  const Double_t AcordePlane = 850.; // distance of the central Acorde detectors to the beam line at y =0
  const Double_t roof = 210.5;       // distance from x =0 to end of magnet roof

  Double_t r[3];
  upperTrack->GetXYZ(r);
  Double_t d[3];
  upperTrack->GetDirection(d);
  Double_t x,z;
  z = r[2] + (d[2]/d[1])*(AcordePlane - r[1]);
  x = r[0] + (d[0]/d[1])*(AcordePlane - r[1]);
  
  if (x > roof) {
    x = r[0] + (d[0]/(d[0]+d[1]))*(AcordePlane+roof-r[0]-r[1]);
    z = r[2] + (d[2]/(d[0]+d[1]))*(AcordePlane+roof-r[0]-r[1]);
  }
  if (x < -roof) {
    x = r[0] + (d[0]/(d[1]-d[0]))*(AcordePlane+roof+r[0]-r[1]);	      
    z = r[2] + (d[2]/(d[1]-d[0]))*(AcordePlane+roof+r[0]-r[1]);
  } 

  fModules->Fill(z, x);
 
}



Long64_t AliTPCcalibCosmic::Merge(TCollection *li) {

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
  Double_t MIPvalue = TMath::Min(funcDoubleGaus->GetParameter(1),funcDoubleGaus->GetParameter(4));

  delete funcDoubleGaus;

  return MIPvalue;

}




void AliTPCcalibCosmic::CalculateBetheParams(TH2F */*hist*/, Double_t * /*initialParam*/) {
  //
  // Not implemented yet
  //
  return;

}


void AliTPCcalibCosmic::BinLogX(TH1 *h) {

  // Method for the correct logarithmic binning of histograms

  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *new_bins = new Double_t[bins + 1];
   
  new_bins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete new_bins;
  
}

AliExternalTrackParam *AliTPCcalibCosmic::Invert(AliExternalTrackParam *input)
{
  //
  // Invert paramerameter  - not covariance yet
  //
  AliExternalTrackParam *output = new AliExternalTrackParam(*input);
  Double_t * param = (Double_t*)output->GetParameter();
  param[0]*=-1;
  param[3]*=-1;
  param[4]*=-1;
  //
  return output;
}

AliExternalTrackParam *AliTPCcalibCosmic::MakeTrack(const AliExternalTrackParam *track0, const AliExternalTrackParam *track1){
  //
  // 
  //
  AliExternalTrackParam *par1R= new AliExternalTrackParam(*track1);
  par1R->Rotate(track0->GetAlpha());
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
  par1R->PropagateTo(track0->GetX(),0); // bz shold be set -
  //if (1){
  //  printf("Print param\n");
  //  track1->Print();
  //  par1R->Print();
  //}
  return par1R;
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
    vecXk(ipar,0) = param1[ipar];
    vecZk(ipar,0) = param2[ipar];
    for (Int_t jpar=0; jpar<5; jpar++){
      covXk(ipar,jpar) = covar1[track1.GetIndex(ipar, jpar)];
      measR(ipar,jpar) = covar2[track2.GetIndex(ipar, jpar)];
    }
  }
  //
  //
  //
  //
  matHk(0,0)=1; matHk(1,1)= 1; matHk(2,2)= 1;  
  matHk(3,3)= 1;    matHk(4,4)= 1;           // vector to measurement
  //
  vecYk = vecZk-matHk*vecXk;                 // Innovation or measurement residual
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;      // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;              //  Optimal Kalman gain
  vecXk += matKk*vecYk;                      //  updated vector 
  mat1(0,0)=1; mat1(1,1)=1; mat1(2,2)=1; mat1(3,3)=1; mat1(4,4)=1;
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

void AliTPCcalibCosmic::ProcessTree(TTree * chainTracklet, AliExternalComparison *comp){
  //
  // Process the debug streamer tree
  // Possible to modify selection criteria
  //
  TTreeSRedirector * cstream = new TTreeSRedirector("cosmicdump.root");
  //AliTPCcalibCosmic *cosmic = this;
  //
  AliExternalTrackParam * tr0 = 0;
  AliExternalTrackParam * tr1 = 0;  
  Int_t npoints =0;
  {
    Int_t entries=chainTracklet->GetEntries();
    for (Int_t i=0; i< entries; i++){
      chainTracklet->GetBranch("Tr0.")->SetAddress(&tr0);
      chainTracklet->GetBranch("Tr1.")->SetAddress(&tr1);
      chainTracklet->GetEntry(i);
      if (!tr0) continue;
      if (!tr1) continue;
      if (tr0->GetY()==0) continue;
      if (tr1->GetY()==0) continue;
      // make a local copy
      AliExternalTrackParam par0(*tr0);
      AliExternalTrackParam par1(*tr1);
      AliExternalTrackParam par1R(*tr1);
      par1R.Rotate(par1.GetAlpha()+TMath::Pi());
      AliExternalTrackParam *par1T = MakeTrack(tr0,tr1);
      if (0) {
	printf("%d\t%d\t\n",i, npoints);
	par1R.Print();
	par1T->Print();
      }
      AliExternalTrackParam par0U=par0;
      AliExternalTrackParam par1U=*par1T;
      //
      UpdateTrack(par0U,*par1T);
      UpdateTrack(par1U,par0);
      //
      //
      if (i%100==0) printf("%d\t%d\tt\n",i, npoints);
      Bool_t accept =   comp->AcceptPair(&par0,par1T);  

      if (1||fStreamLevel>0){
	(*cstream)<<"Tracklet"<<
	  "accept="<<accept<<
	  "tr0.="<<&par0<<       //original track  up
	  "tr1.="<<&par1<<       //original track  down
	  "tr1R.="<<&par1R<<     //track1 rotated to  0 frame 
	  "tr1T.="<<par1T<<      //track1 transformed to the track 0 frame 
	  //
	  "tr0U.="<<&par0U<<     //track 0 updated with track 1
	  "tr1U.="<<&par1U<<     //track 1 updated with track 0 
	  "\n";
      }
      //
      if (accept) {
	npoints++;
	if (comp) comp->Process(&par0,par1T);
      }
      delete par1T;
    }
  }
  delete cstream;
}







/*

void Init(){
  
.x ~/UliStyle.C
.x ~/rootlogon.C
gSystem->Load("libSTAT.so");
gSystem->Load("libANALYSIS");
gSystem->Load("libTPCcalib");
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");

gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
AliXRDPROOFtoolkit tool; 
TChain * chainCosmic = tool.MakeChain("cosmic.txt","Track0",0,1000000);
chainCosmic->Lookup();

TCut cutT("cutT","abs(Tr1.fP[3]+Tr0.fP[3])<0.01");  // OK
TCut cutD("cutD","abs(Tr0.fP[0]+Tr1.fP[0])<2");     // OK
TCut cutP1("cutP1","abs(Tr0.fP[1]-Tr1.fP[1])<20");   // OK
TCut cutPt("cutPt","abs(Tr1.fP[4]+Tr0.fP[4])<0.1&&abs(Tr0.fP[4])+abs(Tr1.fP[4])<10");
TCut cutN("cutN","min(Orig0.fTPCncls,Orig1.fTPCncls)>50");
TCut cutM("cutM","abs(mag)>0.01");
TCut cutA=cutT+cutD+cutPt+cutN+cutP1+"trigger!=16";

TCut cuthpt("abs(Tr0.fP[4])+abs(Tr1.fP[4])<0.2");
TCut cutS("cutS","Orig0.fIp.fP[1]*Orig1.fIp.fP[1]>0");

//
chainCosmic->Draw(">>listELP",cutA,"entryList");
//TEntryList *elist = (TEntryList*)gDirectory->Get("listEL");
//TEntryList *elist = (TEntryList*)gProof->GetOutputList()->At(1);
chainCosmic->SetEntryList(elist);
//
chainCosmic->Draw(">>listV40Z100","abs(d0)<40&&abs(v01)<100","entryList");
TEntryList *elistV40Z100 = (TEntryList*)gDirectory->Get("listV40Z100");
chainCosmic->SetEntryList(elistV40Z100);

//
// Aliases
//
chainCosmic->SetAlias("side","(-1+(Tr0.fP[1]>0)*2)");
chainCosmic->SetAlias("hpt","abs(Tr0.fP[4])<0.2");
chainCosmic->SetAlias("signy","(-1+(Tr0.fP[0]>0)*2)");

chainCosmic->SetAlias("dy","Tr0.fP[0]+Tr1.fP[0]");
chainCosmic->SetAlias("dz","Tr0.fP[1]-Tr1.fP[1]");
chainCosmic->SetAlias("d1pt","Tr0.fP[4]+Tr1.fP[4]");
chainCosmic->SetAlias("dtheta","(Tr0.fP[3]+Tr1.fP[3])");
chainCosmic->SetAlias("dphi","(Tr0.fAlpha-Tr1.fAlpha-pi)");

chainCosmic->SetAlias("mtheta","(Tr0.fP[3]-Tr1.fP[3])*0.5")
chainCosmic->SetAlias("sa","(sin(Tr0.fAlpha+0.))");
chainCosmic->SetAlias("ca","(cos(Tr0.fAlpha+0.))");



chainCosmic->Draw("dy:sqrt(abs(Tr0.fP[4]))>>hisdyA(5,0,1,50,-1,1)","!crossO&&!crossI&&abs(d0)<40&&Tr0.fP[1]>0");
hisdyA->FitSlicesY();
hisdyA_2->SetXTitle("#sqrt{1/p_{t}}");
hisdyA_2->SetYTitle("#sigma_{y}(cm)");
hisdyA_2->SetTitle("Cosmic - Y matching");
hisdyA_2->SetMaximum(0.5);


chainCosmic->Draw("dy:sqrt(abs(Tr0.fP[4]))>>hisdyC(5,0,1,50,-1,1)","!crossO&&!crossI&&abs(d0)<40&&Tr0.fP[1]<0");
hisdyC->FitSlicesY();
hisdyC_2->SetXTitle("#sqrt{1/p_{t}}");
hisdyC_2->SetYTitle("#sigma_{y}(cm)");
hisdyC_2->SetTitle("Cosmic - Y matching");
hisdyC_2->SetMaximum(1);
hisdyC_2->SetMinimum(0);
hisdyC_2->SetMarkerStyle(22);
hisdyA_2->SetMarkerStyle(21);
hisdyC_2->SetMarkerSize(1.5);
hisdyA_2->SetMarkerSize(1.5);
hisdyC_2->Draw();
hisdyA_2->Draw("same");
gPad->SaveAs("~/Calibration/Cosmic/pic/ymatching.gif")

chainCosmic->Draw("dz:sqrt(abs(Tr0.fP[4]))>>hisdzA(5,0,1,50,-1,1)","!crossO&&!crossI&&abs(d0)<40&&Tr0.fP[1]>0");
hisdzA->FitSlicesY();
hisdzA_2->SetXTitle("#sqrt{1/p_{t}}");
hisdzA_2->SetYTitle("#sigma_{z}(cm)");
hisdzA_2->SetTitle("Cosmic - Z matching - A side ");
hisdzA_2->SetMaximum(0.5);

chainCosmic->Draw("dz:sqrt(abs(Tr0.fP[4]))>>hisdzC(5,0,1,50,-1,1)","!crossO&&!crossI&&abs(d0)<40&&Tr0.fP[1]<0");
hisdzC->FitSlicesY();
hisdzC_2->SetXTitle("#sqrt{1/p_{t}}");
hisdzC_2->SetYTitle("#sigma_{z}(cm)");
hisdzC_2->SetTitle("Cosmic - Z matching");
hisdzC_2->SetMaximum(0.5);
hisdzC_2->SetMarkerStyle(22);
hisdzA_2->SetMarkerStyle(21);
hisdzC_2->SetMarkerSize(1.5);
hisdzA_2->SetMarkerSize(1.5);

hisdzC_2->Draw();
hisdzA_2->Draw("same");


//
// PICTURE 1/pt
//
chainCosmic->Draw("d1pt:sqrt(abs(Tr0.fP[4]))>>hisd1ptA(5,0,1,30,-0.1,0.1)","!crossO&&!crossI&&abs(d0)<40&&Tr0.fP[1]>0"+cutM);
hisd1ptA->FitSlicesY();
hisd1ptA_2->SetXTitle("#sqrt{1/p_{t}}");
hisd1ptA_2->SetYTitle("#sigma_{z}(cm)");
hisd1ptA_2->SetTitle("Cosmic - Z matching - A side ");
hisd1ptA_2->SetMaximum(0.5);

chainCosmic->Draw("d1pt:sqrt(abs(Tr0.fP[4]))>>hisd1ptC(5,0,1,30,-0.1,0.1)","!crossO&&!crossI&&abs(d0)<40&&Tr0.fP[1]<0"+cutM);
hisd1ptC->FitSlicesY();
hisd1ptC_2->SetXTitle("#sqrt{1/p_{t}}");
hisd1ptC_2->SetYTitle("#sigma_{1/pt}(1/GeV)");
hisd1ptC_2->SetTitle("Cosmic - 1/pt matching");
hisd1ptC_2->SetMaximum(0.05);
hisd1ptC_2->SetMarkerStyle(22);
hisd1ptA_2->SetMarkerStyle(21);
hisd1ptC_2->SetMarkerSize(1.5);
hisd1ptA_2->SetMarkerSize(1.5);

hisd1ptC_2->Draw();
hisd1ptA_2->Draw("same");
gPad->SaveAs("~/Calibration/Cosmic/pic/1ptmatching.gif")

//
// Theta
//
chainCosmic->Draw("dtheta:sqrt(abs(Tr0.fP[4]))>>hisdthetaA(5,0,1,30,-0.01,0.01)","!crossO&&!crossI&&abs(d0)<40&&Tr0.fP[1]>0");
hisdthetaA->FitSlicesY();
hisdthetaA_2->SetXTitle("#sqrt{1/p_{t}}");
hisdthetaA_2->SetYTitle("#sigma_{#theta}(cm)");
hisdthetaA_2->SetTitle("Cosmic - Z matching - A side ");
hisdthetaA_2->SetMaximum(0.5);

chainCosmic->Draw("dtheta:sqrt(abs(Tr0.fP[4]))>>hisdthetaC(5,0,1,30,-0.01,0.01)","!crossO&&!crossI&&abs(d0)<40&&Tr0.fP[1]<0");
hisdthetaC->FitSlicesY();
hisdthetaC_2->SetXTitle("#sqrt{1/p_{t}}");
hisdthetaC_2->SetYTitle("#sigma_{#theta}(rad)");
hisdthetaC_2->SetTitle("Cosmic - Theta matching");
hisdthetaC_2->SetMaximum(0.01);
hisdthetaC_2->SetMinimum(0.0);
hisdthetaC_2->SetMarkerStyle(22);
hisdthetaA_2->SetMarkerStyle(21);
hisdthetaC_2->SetMarkerSize(1.5);
hisdthetaA_2->SetMarkerSize(1.5);

hisdthetaC_2->Draw();
hisdthetaA_2->Draw("same");
gPad->SaveAs("~/Calibration/Cosmic/pic/thetamatching.gif")
//
// Phi
//
chainCosmic->Draw("dphi:sqrt(abs(Tr0.fP[4]))>>hisdphiA(5,0,1,30,-0.01,0.01)","!crossO&&!crossI&&abs(d0)<40&&Tr0.fP[1]>0");
hisdphiA->FitSlicesY();
hisdphiA_2->SetXTitle("#sqrt{1/p_{t}}");
hisdphiA_2->SetYTitle("#sigma_{#phi}(rad)");
hisdphiA_2->SetTitle("Cosmic - Z matching - A side ");
hisdphiA_2->SetMaximum(0.5);

chainCosmic->Draw("dphi:sqrt(abs(Tr0.fP[4]))>>hisdphiC(5,0,1,30,-0.01,0.01)","!crossO&&!crossI&&abs(d0)<40&&Tr0.fP[1]<0");
hisdphiC->FitSlicesY();
hisdphiC_2->SetXTitle("#sqrt{1/p_{t}}");
hisdphiC_2->SetYTitle("#sigma_{#phi}(rad)");
hisdphiC_2->SetTitle("Cosmic - Phi matching");
hisdphiC_2->SetMaximum(0.01);
hisdphiC_2->SetMinimum(0.0);
hisdphiC_2->SetMarkerStyle(22);
hisdphiA_2->SetMarkerStyle(21);
hisdphiC_2->SetMarkerSize(1.5);
hisdphiA_2->SetMarkerSize(1.5);

hisdphiC_2->Draw();
hisdphiA_2->Draw("same");
gPad->SaveAs("~/Calibration/Cosmic/pic/phimatching.gif")



}


*/


/*
void MatchTheta(){

TStatToolkit toolkit;
Double_t chi2=0;
Int_t    npoints=0;
TVectorD fitParamA0;
TVectorD fitParamA1;
TVectorD fitParamC0;
TVectorD fitParamC1;
TMatrixD covMatrix;


TString fstring="";
// 
fstring+="mtheta++";
fstring+="ca++";
fstring+="sa++";
fstring+="ca*mtheta++";
fstring+="sa*mtheta++";
//
fstring+="side++";
fstring+="side*mtheta++";
fstring+="side*ca++";
fstring+="side*sa++";
fstring+="side*ca*mtheta++";
fstring+="side*sa*mtheta++";


TString *strTheta0 = toolkit.FitPlane(chain,"dtheta",fstring->Data(), "hpt&&!crossI&&!crossO", chi2,npoints,fitParamA0,covMatrix,0.8);
chainCosmic->SetAlias("dtheta0",strTheta0.Data());
strTheta0->Tokenize("+")->Print();


//fstring+="mtheta++";
//fstring+="mtheta^2++";
//fstring+="ca*mtheta^2++";
//fstring+="sa*mtheta^2++";



}

*/




/*
 void PosCorrection()

 

 
 TStatToolkit toolkit;
 Double_t chi2=0;
 Int_t    npoints=0;
 TVectorD fitParam;
 TMatrixD covMatrix;
 
 //Theta
chainCosmic->SetAlias("dthe","(Tr0.fP[3]+Tr1.fP[3])");
chainCosmic->SetAlias("sign","(-1+(Tr0.fP[1]>0)*2)");
chainCosmic->SetAlias("di","(1-abs(Tr0.fP[1])/250)");
//
//
TString strFit="";
//
strFit+="sign++";                              //1
strFit+="Tr0.fP[3]++";                         //2
// 
strFit+="sin(Tr0.fAlpha)*(Tr0.fP[1]>0)++";     //3
strFit+="sin(Tr0.fAlpha)*(Tr0.fP[1]<0)++";     //4
//
strFit+="cos(Tr0.fAlpha)*(Tr0.fP[1]>0)++";     //5
strFit+="cos(Tr0.fAlpha)*(Tr0.fP[1]<0)++";     //6  
//
strFit+="sin(Tr0.fAlpha)*Tr0.fP[3]++";         //7
strFit+="cos(Tr0.fAlpha)*Tr0.fP[3]++";         //8
 
 
 //					    
 TString * thetaParam = toolkit.FitPlane(chain,"dthe", strFit.Data(),"1", chi2,npoints,fitParam,covMatrix,0.8,0,10000)
 
 chainCosmic->SetAlias("corTheta",thetaParam->Data());
 chainCosmic->Draw("dthe:Tr0.fP[1]","","",50000);



*/



/*

void AliTPCcalibCosmic::dEdxCorrection(){
  TCut cutT("cutT","abs(Tr1.fP[3]+Tr0.fP[3])<0.01");  // OK
  TCut cutD("cutD","abs(Tr0.fP[0]+Tr1.fP[0])<2");     // OK
  TCut cutPt("cutPt","abs(Tr1.fP[4]+Tr0.fP[4])<0.1&&abs(Tr0.fP[4])+abs(Tr1.fP[4])<10");
  TCut cutN("cutN","min(Orig0.fTPCncls,Orig1.fTPCncls)>100");
  TCut cutS("cutS","Orig0.fIp.fP[1]*Orig1.fIp.fP[1]>0");
  TCut cutA=cutT+cutD+cutPt+cutN+cutS;


 .x ~/UliStyle.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
 AliXRDPROOFtoolkit tool; 
  TChain * chainCosmic = tool.MakeChain("cosmic.txt","Track0",0,1000000);
  chainCosmic->Lookup();
  
  chainCosmic->Draw(">>listEL",cutA,"entryList");
  TEntryList *elist = (TEntryList*)gDirectory->Get("listEL");
  chainCosmic->SetEntryList(elist);

  .x ~/rootlogon.C
   gSystem->Load("libSTAT.so");
   TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParam;
  TMatrixD covMatrix;
  
  chainCosmic->Draw("Tr0.fP[4]+Tr1.fP[4]","OK"+cutA);
  
  TString strFit;
  strFit+="(Tr0.fP[1]/250)++";
  strFit+="(Tr0.fP[1]/250)^2++";
  strFit+="(Tr0.fP[3])++";
  strFit+="(Tr0.fP[3])^2++";

  TString * ptParam = TStatToolkit::FitPlane(chain,"Tr0.fP[4]+Tr1.fP[4]", strFit.Data(),"1", chi2,npoints,fitParam,covMatrix) 



*/


/*
gSystem->Load("libANALYSIS");
gSystem->Load("libSTAT");
gSystem->Load("libTPCcalib");

TStatToolkit toolkit;
Double_t chi2;
TVectorD fitParam;
TMatrixD covMatrix;
Int_t npoints;
//
TCut cutT("cutT","abs(Tr1.fP[3]+Tr0.fP[3])<0.03");  // OK
TCut cutD("cutD","abs(Tr0.fP[0]+Tr1.fP[0])<5");     // OK
TCut cutPt("cutPt","abs(Tr1.fP[4]+Tr0.fP[4])<0.2&&abs(Tr0.fP[4])+abs(Tr1.fP[4])<10");
TCut cutN("cutN","min(Orig0.fTPCncls,Orig1.fTPCncls)>110");
TCut cutA=cutT+cutD+cutPt+cutN;



TTree * chainCosmic = Track0;


chainCosmic->SetAlias("norm","signalTot0.fElements[3]/signalTot1.fElements[3]");
//
chainCosmic->SetAlias("dr1","(signalTot0.fElements[1]/signalTot0.fElements[3])");
chainCosmic->SetAlias("dr2","(signalTot0.fElements[2]/signalTot0.fElements[3])");
chainCosmic->SetAlias("dr4","(signalTot0.fElements[4]/signalTot0.fElements[3])");
chainCosmic->SetAlias("dr5","(signalTot0.fElements[5]/signalTot0.fElements[3])");

TString fstring="";  
fstring+="dr1++";
fstring+="dr2++";
fstring+="dr4++";
fstring+="dr5++";
//
fstring+="dr1*dr2++";
fstring+="dr1*dr4++";
fstring+="dr1*dr5++";
fstring+="dr2*dr4++";
fstring+="dr2*dr5++";
fstring+="dr4*dr5++";



TString *strqdedx = toolkit.FitPlane(chain,"norm",fstring->Data(), cutA, chi2,npoints,fitParam,covMatrix,-1,0,200000);
  
chainCosmic->SetAlias("corQT",strqdedx->Data());

*/


/*
  chainCosmic->SetProof(kTRUE);
  chainCosmic->Draw("Seed0.CookdEdxNorm(0,0.6,1,0,159,0,kTRUE,kTRUE):Seed0.CookdEdxNorm(0,0.6,1,0,159,0,kFALSE,kTRUE)",""+cutA,"",100000);


chainCosmic->Draw("Seed0.CookdEdxNorm(0,0.6,1,0,159,0,kTRUE,kTRUE)/Seed1.CookdEdxNorm(0,0.6,1,0,159,0,kTRUE,kTRUE)>>his(100,0.5,1.5)","min(Orig0.fTPCncls,Orig1.fTPCncls)>130"+cutA,"",50000);

*/


/*
chainCosmic->Draw("Tr0.fP[1]-Tr1.fP[1]:sqrt(abs(Tr0.fP[4]))>>hisdzA(5,0,1,50,-1,1)","!crossO&&!crossI&&abs(d0)<40&&Tr0.fP[1]>0&&abs(mag)>0.1"+cutA); 

TGraph *grdzA = (TGraph*)gProof->GetOutputList()->At(1)->Clone();



 
*/




