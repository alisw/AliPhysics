#include "AliHBTTrackPoints.h"
#include "AliTPCtrack.h"
#include <TMath.h>
#include "AliTrackReference.h"
#include <TClonesArray.h>

ClassImp(AliHBTTrackPoints)

AliHBTTrackPoints::AliHBTTrackPoints():
 fN(0),
 fX(0x0),
 fY(0x0),
 fZ(0x0)
{
  //constructor
}

AliHBTTrackPoints::AliHBTTrackPoints(Int_t n, AliTPCtrack* track, Float_t dr, Float_t r0):
 fN(n),
 fX(new Float_t[fN]),
 fY(new Float_t[fN]),
 fZ(new Float_t[fN])
{
  //constructor
  //r0 starting radius
  //dr - calculate points every dr cm, default every 30cm
  if (track == 0x0)
   {
     Error("AliHBTTrackPoints","TPC track is null");
     fN = 0;
     delete [] fX;
     delete [] fY;
     delete [] fZ;
     fX = fY = fZ = 0x0;
     return;
   }
  track->PropagateTo(r0);

 //* This formation is now fixed in the following way:                         *
 //*      external param0:   local Y-coordinate of a track (cm)                *
 //*      external param1:   local Z-coordinate of a track (cm)                *
 //*      external param2:   local sine of the track momentum azimuth angle    *
 //*      external param3:   tangent of the track momentum dip angle           *
 //*      external param4:   1/pt (1/(GeV/c))                                  *
    
  Double_t xk;
  Double_t par[5];
  track->GetExternalParameters(xk,par);     //get properties of the track
  
  Double_t x = track->GetX();
  Double_t y = track->GetY();
  Double_t z0 = track->GetZ();
      
  Double_t phi0local = TMath::ATan2(y,x);
  Double_t phi0global = phi0local + track->GetAlpha();
  
  if (phi0local<0) phi0local+=2*TMath::Pi();
  if (phi0local>=2.*TMath::Pi()) phi0local-=2*TMath::Pi();
  
  if (phi0global<0) phi0global+=2*TMath::Pi();
  if (phi0global>=2.*TMath::Pi()) phi0global-=2*TMath::Pi();
  
  Double_t r = TMath::Hypot(x,y);
  

  if (GetDebug()) 
    Info("AliHBTTrackPoints","Radius0 %f, Real Radius %f",r0,r); 
  
  if (GetDebug()) 
    Info("AliHBTTrackPoints","Phi Global at first padraw %f, Phi locat %f",phi0global,phi0local);
  Double_t c=track->GetC();
  Double_t eta = track->GetEta();
  
  //this calculattions are assuming new (current) model 
  Double_t tmp = c*x - eta;
  tmp = 1. - tmp*tmp;
  tmp = c*y + TMath::Sqrt(tmp);
  Double_t dca=(TMath::Hypot(eta,tmp) - 1. )/TMath::Abs(c);

  //Here is old model Cold=Cnew/2.
  Double_t dcasq = dca*dca;
  Double_t c2 = c/2.;
  Double_t cst1 = (1.+c2*dca)*dca;//first constant
  Double_t cst2 = 1. + 2.*c2*dca;//second constant
      
  Double_t factorPhi0 = TMath::ASin((c2*r + cst1/r)/cst2);
  Double_t factorZ0 = TMath::ASin(c2*TMath::Sqrt((r*r-dcasq)/cst2))*par[3]/c2;
 
  for(Int_t i = 0; i<fN; i++)
   {
     Double_t rc = r0 + i*dr;
     Double_t factorPhi = TMath::ASin( (c2*rc + cst1/rc)/cst2 );//factor phi od rc
     Double_t phi = phi0global + factorPhi - factorPhi0;
     
     Double_t factorZ = TMath::ASin(c2*TMath::Sqrt((rc*rc-dcasq)/cst2))*par[3]/c2;
     fZ[i] = z0 + factorZ - factorZ0;
     fX[i] = rc*TMath::Cos(phi);
     fY[i] = rc*TMath::Sin(phi);
     
     if ( GetDebug() )
      {
        Info("AliHBTTrackPoints","X %f Y %f Z %f R asked %f R obtained %f",
             fX[i],fY[i],fZ[i],rc,TMath::Hypot(fX[i],fY[i]));
      }
   }
}


AliHBTTrackPoints::~AliHBTTrackPoints()
{
  //destructor
  delete [] fX;
  delete [] fY;
  delete [] fZ;
}

void AliHBTTrackPoints::PositionAt(Int_t n, Float_t &x,Float_t &y,Float_t &z)
{
  //returns position at point n
  if ((n<0) || (n>fN))
   {
     Error("PositionAt","Point %d out of range",n);
     return;
   }
  
  x = fX[n];
  y = fY[n];
  z = fZ[n];
  if ( GetDebug() )
    {
      Info("AliHBTTrackPoints","n %d; X %f; Y %f; Z %f",n,x,y,z);
    }

}


Double_t AliHBTTrackPoints::AvarageDistance(const AliHBTTrackPoints& tr)
{
  //returns the aritmethic avarage distance between two tracks
  if (fN != tr.fN)
   {
     Error("AvarageDistance","Number of points is not equal");
     return -1;
   }
  if ( (fN <= 0) || (tr.fN <=0) )
   {
     Error("AvarageDistance","One of tracks is empty");
     return -1;
   }
   
  Double_t sum;
  for (Int_t i = 0; i<fN; i++)
   {
      if (fR[i] != tr.fR[i])
       {
         Error("AvarageDistance","Different radii");
         return -1;
       }
      Double_t dx = fX-tr.fX;
      Double_t dy = fY-tr.fY;
      Double_t dz = fZ-tr.fZ;
      sum+=TMath::Sqrt(dx*dx + dy*dy + dz*dz);
   }
  return sum/((Double_t)fN);
}

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliTPCtrack.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "AliMagF.h"



void AliHBTTrackPoints::tp(Int_t entr)
{
  delete gAlice;
  gAlice = 0x0;
  AliRunLoader* rl = AliRunLoader::Open();
  AliLoader* l = rl->GetLoader("TPCLoader");
  rl->LoadgAlice();
  AliKalmanTrack::SetConvConst(100/0.299792458/0.2/rl->GetAliRun()->Field()->Factor());
  l->LoadTracks();
  AliTPCtrack* t = new AliTPCtrack();
  TBranch* b=l->TreeT()->GetBranch("tracks");
  b->SetAddress(&t);
  l->TreeT()->GetEntry(entr);  
  Int_t N = 160;
  AliHBTTrackPoints* tp = new AliHBTTrackPoints(N,t,1.);
  
  TH2D* hxy = new TH2D("hxy","hxy",1000,150,250,1000,150,250);
  TH2D* hxyt = new TH2D("hxyt","hxyt",1000,150,250,1000,150,250);
  TH2D* hxyTR = new TH2D("hxyTR","hxyTR",1000,150,250,1000,150,250);

  TH2D* hxz = new TH2D("hxz","hxz",1000,150,250,1000,151,251);
  TH2D* hxzt = new TH2D("hxzt","hxzt",1000,150,250,1000,151,251);
  TH2D* hxzTR = new TH2D("hxzTR","hxzTR",1000,150,250,1000,151,251);

  hxyt->SetDirectory(0x0);   
  hxy->SetDirectory(0x0);   
  hxyTR->SetDirectory(0x0);   
  
  hxzt->SetDirectory(0x0);   
  hxz->SetDirectory(0x0);   
  hxzTR->SetDirectory(0x0);   

  Float_t x,y,z;
   
  for (Int_t i = 0;i<N;i++)
   {  
      Double_t r = 84.1+i;
      tp->PositionAt(i,x,y,z);
      hxy->Fill(x,y);
      hxz->Fill(x,z);
      printf("Rdemanded %f\n",r);
      printf("tpx %f tpy %f tpz %f Rt =%f\n", x,y,z,TMath::Hypot(x,y));
      
      //BUT they are local!!!!
      t->PropagateTo(r);
//      Double_t phi = t->Phi();
      Double_t rl = TMath::Hypot(t->GetX(),t->GetY());//real radius 
      
      Double_t alpha = t->GetAlpha();
      Double_t salpha = TMath::Sin(alpha);
      Double_t calpha = TMath::Cos(alpha);
      x = t->GetX()*calpha - t->GetY()*salpha;
      y = t->GetX()*salpha + t->GetY()*calpha;
      z = t->GetZ();
      
      printf("tx %f ty %f tz %f Rt = %f R from XY %f\n",x,y,z,TMath::Hypot(x,y),rl);
      
      printf("tpz - tz %f\n",z-t->GetZ());
      printf("\n");
      hxyt->Fill(x,y);
      hxzt->Fill(x,z);
      
   }

  rl->LoadTrackRefs();
  TTree* treeTR = rl->TreeTR();
  b = treeTR->GetBranch("TPC");
  
  TClonesArray* trackrefs = new TClonesArray("AliTrackReference", 100);
  AliTrackReference* tref;
  b->SetAddress(&trackrefs);
  
  Int_t tlab = TMath::Abs(t->GetLabel());

  Int_t netr = (Int_t)treeTR->GetEntries();
  printf("Found %d entries in TR tree\n",netr);
  
  for (Int_t e = 0; e < netr; e++)
   { 
     treeTR->GetEntry(e);
     tref = (AliTrackReference*)trackrefs->At(0);
     if (tref == 0x0) continue;
     if (tref->GetTrack() != tlab) continue;
    
     printf("Found %d entries in TR array\n",trackrefs->GetEntries());
     
     for (Int_t i = 0; i < trackrefs->GetEntries(); i++)
      {
        tref = (AliTrackReference*)trackrefs->At(i);
        if (tref->GetTrack() != tlab) continue;
        x = tref->X();
        y = tref->Y();
        z = tref->Z();
        printf("Track Ref: x %f y %f z %f\n",tref->X(),tref->Y(),tref->Z());

        hxzTR->Fill(x,z);
        hxyTR->Fill(x,y);
        for (Int_t j = 1; j < 10; j++)
         {
            hxyTR->Fill(x, y+j*0.1);
            hxyTR->Fill(x, y-j*0.1);
            hxyTR->Fill(x+j*0.1,y);
            hxyTR->Fill(x-j*0.1,y);

            hxzTR->Fill(x,z-j*0.1);
            hxzTR->Fill(x,z+j*0.1);
            hxzTR->Fill(x-j*0.1,z);
            hxzTR->Fill(x+j*0.1,z);
         }
      }
      break; 
    }  
  hxz->Draw("");
//  hxzt->Draw("same");
  hxzTR->Draw("same");
  
  delete rl;
}
