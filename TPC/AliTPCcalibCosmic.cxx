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
  fDeDx(0)
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
  fDeDx(0)
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
  AliInfo("Non Defualt Constructor");  
}

AliTPCcalibCosmic::~AliTPCcalibCosmic(){
  //
  //
  //
}


void AliTPCcalibCosmic::Process(AliESDEvent *event) {
    
  if (!event) {
    Printf("ERROR: ESD not available");
    return;
  }
  
  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!ESDfriend) {
   Printf("ERROR: ESDfriend not available");
   return;
  }

  printf("Hallo world: Im here\n");
  
  Int_t n=event->GetNumberOfTracks(); 
  fHistNTracks->Fill(n);
  
  //track loop
  for (Int_t i=0;i<n;++i) { 
   AliESDtrack *track = event->GetTrack(i); 
   fClusters->Fill(track->GetTPCNcls());
   
   AliExternalTrackParam * trackIn = new AliExternalTrackParam(*track->GetInnerParam());
   
   AliESDfriendTrack *friendTrack = ESDfriend->GetTrack(i);
   TObject *calibObject;
   AliTPCseed *seed = 0;
   for (Int_t l=0;calibObject=friendTrack->GetCalibObject(l);++l) {
    if (seed=dynamic_cast<AliTPCseed*>(calibObject)) break;
   }
 
   if (seed && track->GetTPCNcls() > 80) fDeDx->Fill(trackIn->GetP(), seed->CookdEdxNorm(0.05,0.45,0));

    
  }
  
  // dE/dx,pt and ACORDE study --> studies which need the pair selection  
  if (n > 2) return;
  
  for (Int_t i=0;i<n;++i) {
    AliESDtrack *track1 = event->GetTrack(i);
     
    Double_t d1[3];
    track1->GetDirection(d1);
    
    for (Int_t j=i+1;j<n;++j) {
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
      
      printf("My stream level=%d\n",fStreamLevel);

      if (fStreamLevel>0){
	TTreeSRedirector * cstream =  GetDebugStreamer();
	printf("My stream=%p\n",cstream);
	if (cstream) {
	  (*cstream) << "Track" <<
	    "Track1.=" << track1 <<      // original track 1
	    "Track2.=" << track2 <<      // original track2
	    "\n";
	}
      }
  //     AliExternalTrackParam * trackOut = new AliExternalTrackParam(*track2->GetOuterParam());
//       AliTracker::PropagateTrackTo(trackOut,850.,105.658,30);
//       delete trackOut;
      


      

      break;            
     }     
    }
  }
  
  
  
  
}    


Long64_t AliTPCcalibCosmic::Merge(TCollection *li) {

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

