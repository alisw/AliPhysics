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

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"

#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"

#include "AliTracker.h"
#include "AliMagFMaps.h"

#include "AliLog.h"

#include "AliTPCcalibCosmic.h"

#include "TTreeStream.h"
#include "AliTPCTracklet.h"

ClassImp(AliTPCcalibCosmic)


AliTPCcalibCosmic::AliTPCcalibCosmic() 
  :AliTPCcalibBase(),
   fHistNTracks(0),
   fClusters(0),
   fModules(0),
   fHistPt(0),
   fPtResolution(0),
   fDeDx(0),
   fCutMaxD(5),        // maximal distance in rfi ditection
   fCutTheta(0.03),    // maximal distan theta
   fCutMinDir(-0.99)   // direction vector products
{  
  AliInfo("Defualt Constructor");  
}


AliTPCcalibCosmic::AliTPCcalibCosmic(const Text_t *name, const Text_t *title) 
  :AliTPCcalibBase(),
   fHistNTracks(0),
   fClusters(0),
   fModules(0),
   fHistPt(0),
   fPtResolution(0),
   fDeDx(0),
   fCutMaxD(5),        // maximal distance in rfi ditection
   fCutTheta(0.03),    // maximal distan theta
   fCutMinDir(-0.99)   // direction vector products
{  
  SetName(name);
  SetTitle(title);
  AliMagFMaps * field = new AliMagFMaps("dummy1", "dummy2",0,5,0);
  AliTracker::SetFieldMap(field, kTRUE);  
  fHistNTracks = new TH1F("ntracks","Number of Tracks per Event",501,-0.5,500.5);
  fClusters = new TH1F("signal","Number of Clusters per track",160,0,160);
  fModules = new TH2F("sector","Acorde hits; z (cm); x(cm)",1200,-1200,1200,600,-1000,1000);
  fHistPt = new TH1F("Pt","Pt distribution",2000,0,50);  
  fPtResolution = new TH1F("PtResolution","Pt resolution",100,-50,50);
  fDeDx = new TH2F("DeDx","dEdx",500,0.01,20.,500,0.,500);
  BinLogX(fDeDx);
  AliInfo("Non Default Constructor");  
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
  FindPairs(event);

  if (GetDebugLevel()>1) printf("Hallo world: Im here\n");
  Int_t ntracks=event->GetNumberOfTracks(); 
  fHistNTracks->Fill(ntracks);
  TObjArray  tpcSeeds(ntracks);
  if (ntracks==0) return;
  //
  //track loop
  //
  for (Int_t i=0;i<ntracks;++i) { 
   AliESDtrack *track = event->GetTrack(i); 
   fClusters->Fill(track->GetTPCNcls());   
   AliExternalTrackParam * trackIn = new AliExternalTrackParam(*track->GetInnerParam());
   
   AliESDfriendTrack *friendTrack = ESDfriend->GetTrack(i);
   TObject *calibObject;
   AliTPCseed *seed = 0;
   for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
     if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
   }
   if (seed) tpcSeeds.AddAt(seed,i);
   if (seed && track->GetTPCNcls() > 80) fDeDx->Fill(trackIn->GetP(), seed->CookdEdxNorm(0.05,0.45,0)); 
  }
  if (ntracks<2) return;


  // dE/dx,pt and ACORDE study --> studies which need the pair selection    
  for (Int_t i=0;i<ntracks;++i) {
    AliESDtrack *track1 = event->GetTrack(i);
     
    Double_t d1[3];
    track1->GetDirection(d1);
    
    for (Int_t j=i+1;j<ntracks;++j) {
     AliESDtrack *track2 = event->GetTrack(j);   
     Double_t d2[3];
     track2->GetDirection(d2);
       
     if (d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2] < -0.999) {
     
      /*___________________________________ Pt resolution ________________________________________*/
      if (track1->Pt() != 0 && track1->GetTPCNcls() > 80 && track2->GetTPCNcls() > 80) {
       Double_t res = (track1->Pt() - track2->Pt());
       res = res/(2*(track1->Pt() + track2->Pt()));
       fPtResolution->Fill(100*res);
      }
      
      /*_______________________________ Propagation to ACORDE ___________________________________*/
      const Double_t AcordePlane = 850.; //distance of the central Acorde detectors to the beam line at y =0
      const Double_t roof = 210.5; // distance from x =0 to end of magnet roof
     
      if (d1[1] > 0 && d2[1] < 0 && track1->GetTPCNcls() > 50) {        
       Double_t r[3];
       track1->GetXYZ(r);
       Double_t x,z;
       z = r[2] + (d1[2]/d1[1])*(AcordePlane - r[1]);
       x = r[0] + (d1[0]/d1[1])*(AcordePlane - r[1]);
       
       if (x > roof) {
        x = x - (x-roof)/(1 + TMath::Abs(TMath::Tan(track1->Phi())));
        z = z - TMath::Abs(TMath::Tan(track1->Phi()))/TMath::Abs(TMath::Tan(track1->Theta()))*(x-roof)/(1 + TMath::Abs(TMath::Tan(track1->Phi())));
       }
       if (x < -roof) {
        x = x - (x+roof)/(1 + TMath::Abs(TMath::Tan(track1->Phi())));
        z = z -  TMath::Abs(TMath::Tan(track1->Phi()))/TMath::Abs(TMath::Tan(track1->Theta()))*(x+roof)/(1 + TMath::Abs(TMath::Tan(track1->Phi())));
       }
       
       fModules->Fill(z, x);
      }
      
      if (d2[1] > 0 && d1[1] < 0 && track2->GetTPCNcls() > 50) {
       Double_t r[3];
       track2->GetXYZ(r);
       Double_t x,z;
       z = r[2] + (d2[2]/d2[1])*(AcordePlane - r[1]);
       x = r[0] + (d2[0]/d2[1])*(AcordePlane - r[1]);
       
       if (x > roof) {
        x = x - (x-roof)/(1 + TMath::Abs(TMath::Tan(track2->Phi())));
        z = z - TMath::Abs(TMath::Tan(track2->Phi()))/TMath::Abs(TMath::Tan(track2->Theta()))*(x-roof)/(1 + TMath::Abs(TMath::Tan(track2->Phi())));  
       }
       if (x < -roof) {
        x = x - (x+roof)/(1 + TMath::Abs(TMath::Tan(track2->Phi())));
	z = z -  TMath::Abs(TMath::Tan(track2->Phi()))/TMath::Abs(TMath::Tan(track2->Theta()))*(x+roof)/(1 + TMath::Abs(TMath::Tan(track2->Phi())));
       }       
       
       fModules->Fill(z, x);
      }
      
  //     AliExternalTrackParam * trackOut = new AliExternalTrackParam(*track2->GetOuterParam());
//       AliTracker::PropagateTrackTo(trackOut,850.,105.658,30);
//       delete trackOut;
      


      

      break;            
     }     
    }
  }
  
  
  
  
}    

void AliTPCcalibCosmic::FindPairs(AliESDEvent *event) {
  //
  // Find cosmic pairs
  //
  //
  if (GetDebugLevel()>1) printf("Hallo world: Im here\n");
  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  Int_t ntracks=event->GetNumberOfTracks(); 
  fHistNTracks->Fill(ntracks);
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
   AliExternalTrackParam * trackIn = new AliExternalTrackParam(*track->GetInnerParam());
   
   AliESDfriendTrack *friendTrack = ESDfriend->GetTrack(i);
   TObject *calibObject;
   AliTPCseed *seed = 0;
   for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
     if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
   }
   if (seed) tpcSeeds.AddAt(seed,i);
   if (seed && track->GetTPCNcls() > 80) fDeDx->Fill(trackIn->GetP(), seed->CookdEdxNorm(0.05,0.45,0)); 
  }
  if (ntracks<2) return;
  //
  // Find pairs
  //
  for (Int_t i=0;i<ntracks;++i) {
    AliESDtrack *track0 = event->GetTrack(i);     
    Double_t d1[3];
    track0->GetDirection(d1);    
    for (Int_t j=i+1;j<ntracks;++j) {
     AliESDtrack *track1 = event->GetTrack(j);   
     Double_t d2[3];
     track1->GetDirection(d2);
      printf("My stream level=%d\n",fStreamLevel);
      AliTPCseed * seed0 = (AliTPCseed*) tpcSeeds.At(i);
      AliTPCseed * seed1 = (AliTPCseed*) tpcSeeds.At(j);
      if (! seed0) continue;
      if (! seed1) continue;
      Float_t dedx0 = seed0->CookdEdxNorm(0.05,0.55,0);
      Float_t dedx1 = seed1->CookdEdxNorm(0.05,0.55,0);
      Float_t dir = (d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]);
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
      //
      Double_t xyz0[3];//,pxyz0[3];
      Double_t xyz1[3];//,pxyz1[3];
      param0.GetXYZ(xyz0);
      param1.GetXYZ(xyz1);
      Bool_t isPair = IsPair(&param0,&param1);
      //
      if (fStreamLevel>0){
	TTreeSRedirector * cstream =  GetDebugStreamer();
	printf("My stream=%p\n",(void*)cstream);
	if (cstream) {
	  (*cstream) << "Track0" <<
	    "dir="<<dir<<               //  direction
	    "OK="<<isPair<<             // will be accepted
	    "b0="<<b0<<                 //  propagate status
	    "b1="<<b1<<                 //  propagate status
	    "Orig0.=" << track0 <<      //  original track  0
	    "Orig1.=" << track1 <<      //  original track  1
	    "Tr0.="<<&param0<<          //  track propagated to the DCA 0,0
	    "Tr1.="<<&param1<<          //  track propagated to the DCA 0,0	   
	    "v00="<<dvertex0[0]<<       //  distance using kalman
	    "v01="<<dvertex0[1]<<       // 
	    "v10="<<dvertex1[0]<<       //
	    "v11="<<dvertex1[1]<<       // 
	    "d0="<<d0<<                 //  linear distance to 0,0
	    "d1="<<d1<<                 //  linear distance to 0,0
	    //
	    "x00="<<xyz0[0]<<
	    "x01="<<xyz0[1]<<
	    "x02="<<xyz0[2]<<
	    //
	    "x10="<<xyz1[0]<<
	    "x11="<<xyz1[1]<<
	    "x12="<<xyz1[2]<<
	    //
	    "Seed0.=" << track0 <<      //  original seed 0
	    "Seed1.=" << track1 <<      //  original seed 1
	    "dedx0="<<dedx0<<           //  dedx0
	    "dedx1="<<dedx1<<           //  dedx1
	    "\n";
	}
      }      
    }
  }  
}    





Long64_t AliTPCcalibCosmic::Merge(TCollection */*li*/) {
  
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
  if (TMath::Abs(p0[0]+p1[0])>fCutMaxD)  return kFALSE;
  Double_t d0[3], d1[3];
  tr0->GetDirection(d0);    
  tr1->GetDirection(d1);       
  if (d0[0]*d1[0] + d0[1]*d1[1] + d0[2]*d1[2] >fCutMinDir) return kFALSE;
  //
  return kTRUE;  
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

