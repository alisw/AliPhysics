#include "AliHBTTrackPoints.h"
//_________________________________
////////////////////////////////////////////////////////////
//                                                        //
// class AliHBTTrackPoints                                //
//                                                        //
// used by Anti-Merging cut                               //
// contains set of poits the lay on track trajectory      //
// according to reconstructed track parameters -          //
// NOT CLUSTERS POSITIONS!!!                              //
// Anti-Merging cut is applied only on tracks coming from //
// different events (that are use to fill deniminators)   //
//                                                        //
////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TFile.h>
#include <TMath.h>

#include "AliESDtrack.h"
#include "AliTPCtrack.h"
#include "AliTrackReference.h"
#include "AliITStrackV2.h"

ClassImp(AliHBTTrackPoints)

Int_t AliHBTTrackPoints::fgDebug = 0;
AliHBTTrackPoints::AliHBTTrackPoints():
 fN(0),
 fX(0x0),
 fY(0x0),
 fZ(0x0)
{
  //constructor
}
/***************************************************************/

AliHBTTrackPoints::AliHBTTrackPoints(AliHBTTrackPoints::ETypes type, AliESDtrack* track):
 fN(0),
 fX(0x0),
 fY(0x0),
 fZ(0x0)
{
  //constructor
  switch (type)
   {
     case kITS:
       //Used only in non-id analysis
       fN = 6;
       fX = new Float_t[fN];
       fY = new Float_t[fN];
       fZ = new Float_t[fN];
       MakeITSPoints(track);
       break;

     default:
       Info("AliHBTTrackPoints","Not recognized type");
   }
   
}
/***************************************************************/

AliHBTTrackPoints::AliHBTTrackPoints(Int_t n, AliESDtrack* track, Float_t mf, Float_t dr, Float_t r0):
 fN(n),
 fX(new Float_t[fN]),
 fY(new Float_t[fN]),
 fZ(new Float_t[fN])
{
  //constructor
  //mf - magnetic field in kG - needed to calculated curvature out of Pt
  //r0 - starting radius
  //dr - calculate points every dr cm, default every 30cm
  if (track == 0x0)
   {
     Error("AliHBTTrackPoints","ESD track is null");
     fN = 0;
     delete [] fX;
     delete [] fY;
     delete [] fZ;
     fX = fY = fZ = 0x0;
     return;
   }

  if ( ((track->GetStatus() & AliESDtrack::kTPCrefit) == kFALSE)&& 
       ((track->GetStatus() & AliESDtrack::kTPCin) == kFALSE)  )
   {
     //could happend: its stand alone tracking
     if (GetDebug() > 3) 
       Warning("AliHBTTrackPoints","This ESD track does not contain TPC information");
       
     fN = 0;
     delete [] fX;
     delete [] fY;
     delete [] fZ;
     fX = fY = fZ = 0x0;
     
     return;
   }
  
  Double_t x;
  Double_t par[5];
  track->GetInnerExternalParameters(x,par);     //get properties of the track
  if (par[4] == 0)
   {
     Error("AliHBTTrackPoints","This ESD track seem not to contain TPC information (curv is 0)");
     return;
   }

  if (mf == 0.0)
   {
     Error("AliHBTTrackPoints","Zero Magnetic field passed as parameter.");
     return;
   }
   
  Double_t alpha = track->GetInnerAlpha();
  Double_t cc = 1000./0.299792458/mf;//conversion constant
  Double_t c=par[4]/cc;
  
  MakePoints(dr,r0,x,par,c,alpha);
  
}
/***************************************************************/

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
    
  Double_t x = 0;
  Double_t par[5];
  track->GetExternalParameters(x,par);     //get properties of the track
  
  Double_t alpha = track->GetAlpha();
  Double_t c=track->GetC();
  MakePoints(dr,r0,x,par,c,alpha);
}  
/***************************************************************/

AliHBTTrackPoints::~AliHBTTrackPoints()
{
  //destructor
  delete [] fX;
  delete [] fY;
  delete [] fZ;
}
/***************************************************************/

void AliHBTTrackPoints::MakePoints( Float_t dr, Float_t r0, Double_t x, Double_t* par, Double_t c, Double_t alpha)
{
  //Calculates points starting at radius r0
  //spacing every dr (in radial direction)
  // according to track parameters
  // x - position in sector local reference frame. x i parallel to R and sector is symmetric with respect to x
  // par - track external parameters; array with 5 elements;  look at AliTPCtrack.h or AliESDtrack.h for their meaning
  // c - track curvature
  // alpha - sector's rotation angle (phi) == angle needed for local to global transformation
  
  Double_t y = par[0];
  Double_t z0 = par[1];
      
  Double_t phi0local = TMath::ATan2(y,x);
  Double_t phi0global = phi0local + alpha;
  
  if (phi0local<0) phi0local+=2*TMath::Pi();
  if (phi0local>=2.*TMath::Pi()) phi0local-=2*TMath::Pi();
  
  if (phi0global<0) phi0global+=2*TMath::Pi();
  if (phi0global>=2.*TMath::Pi()) phi0global-=2*TMath::Pi();
  
  Double_t r = TMath::Hypot(x,y);
  

  if (GetDebug() > 9) 
    Info("AliHBTTrackPoints","Radius0 %f, Real Radius %f",r0,r); 
  
  if (GetDebug() > 5) 
    Info("AliHBTTrackPoints","Phi Global at first padraw %f, Phi locat %f",phi0global,phi0local);
    
  Double_t eta = x*c - par[2] ;//par[2] = fX*C - eta; eta==fP2 ; C==fP4
  
  //this calculattions are assuming new (current) model 
  Double_t tmp = par[2];
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
     Double_t ftmp = (c2*rc + cst1/rc)/cst2;
     if (ftmp > 1.0) 
      {
        if (GetDebug() > 1) 
          Warning("AliHBTTrackPoints","ASin argument > 1 %f:",ftmp);
        ftmp=1.0;
      }
     else if (ftmp < -1.0) 
      {
        if (GetDebug() > 1) 
          Warning("AliHBTTrackPoints","ASin argument < -1 %f:",ftmp);
        ftmp=-1.0;
      }
      
     Double_t factorPhi = TMath::ASin( ftmp );//factor phi od rc
     Double_t phi = phi0global + factorPhi - factorPhi0;
     
     ftmp = (rc*rc-dcasq)/cst2;
     if (ftmp < 0.0) 
      {
        if (GetDebug() > 1) 
          Warning("AliHBTTrackPoints","Sqrt argument < 0: %f",ftmp);
        ftmp=0.0;
      }
     
     ftmp = c2*TMath::Sqrt(ftmp);
     if (ftmp > 1.0) 
      {
        if (GetDebug() > 1) 
          Warning("AliHBTTrackPoints","ASin argument > 1: %f",ftmp);
        ftmp=1.0;
      }
     else if (ftmp < -1.0) 
      {
        if (GetDebug() > 1) 
          Warning("AliHBTTrackPoints","ASin argument < -1: %f",ftmp);
        ftmp=-1.0;
      }
     Double_t factorZ = TMath::ASin(ftmp)*par[3]/c2;
     fZ[i] = z0 + factorZ - factorZ0;
     fX[i] = rc*TMath::Cos(phi);
     fY[i] = rc*TMath::Sin(phi);
     
     if ( GetDebug() > 2 )
      {
        Info("AliHBTTrackPoints","X %f Y %f Z %f R asked %f R obtained %f",
             fX[i],fY[i],fZ[i],rc,TMath::Hypot(fX[i],fY[i]));
      }
   }
}
/***************************************************************/

void AliHBTTrackPoints::MakeITSPoints(AliESDtrack* track)
{
//Calculates points in ITS
// z=R*Pz/Pt
 AliITStrackV2 itstrack(*track,kTRUE);
 Double_t x,y,z;
 static const Double_t r[6] = {4.0, 7.0, 14.9, 23.8, 39.1, 43.6};
 for (Int_t i = 0; i < 6; i++)
  {
    itstrack.GetGlobalXYZat(r[i],x,y,z);
    fX[i] = x;
    fY[i] = y;
    fZ[i] = z;
//    Info("MakeITSPoints","X %f Y %f Z %f R asked %f R obtained %f",
//             fX[i],fY[i],fZ[i],r[i],TMath::Hypot(fX[i],fY[i]));
  }   
 
}

/***************************************************************/
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
  if ( GetDebug() > 1 )
    {
      Info("AliHBTTrackPoints","n %d; X %f; Y %f; Z %f",n,x,y,z);
    }
}
/***************************************************************/

void AliHBTTrackPoints::Move(Float_t x, Float_t y, Float_t z)
{
//Moves all points about vector
 for (Int_t i = 0; i<fN; i++)
   {
     fX[i]+=x;
     fY[i]+=y;
     fZ[i]+=z;
   }   
}
/***************************************************************/

Double_t AliHBTTrackPoints::AvarageDistance(const AliHBTTrackPoints& tr)
{
  //returns the aritmethic avarage distance between two tracks
//  Info("AvarageDistance","Entered");
  if ( (fN <= 0) || (tr.fN <=0) )
   {
     if (GetDebug()) Warning("AvarageDistance","One of tracks is empty");
     return -1;
   }

  if (fN != tr.fN)
   {
     Warning("AvarageDistance","Number of points is not equal");
     return -1;
   }
   
  Double_t sum = 0;
  for (Int_t i = 0; i<fN; i++)
   {
     if (GetDebug()>9)
      {
//       Float_t r1sq = fX[i]*fX[i]+fY[i]*fY[i];
//       Float_t r2sq = tr.fX[i]*tr.fX[i]+tr.fY[i]*tr.fY[i];
       Float_t r1sq = TMath::Hypot(fX[i],fY[i]);
       Float_t r2sq = TMath::Hypot(tr.fX[i],tr.fY[i]);
       Info("AvarageDistance","radii: %f %f",r1sq,r2sq);
      } 

      
     Double_t dx = fX[i]-tr.fX[i];
     Double_t dy = fY[i]-tr.fY[i];
     Double_t dz = fZ[i]-tr.fZ[i];
     sum+=TMath::Sqrt(dx*dx + dy*dy + dz*dz);
     
     if (GetDebug()>1)
      {
       Info("AvarageDistance","Diff: x ,y z: %f , %f, %f",dx,dy,dz);
       Info("AvarageDistance","xxyyzz %f %f %f %f %f %f",
            fX[i],tr.fX[i],fY[i],tr.fY[i],fZ[i],tr.fZ[i]);
      } 
   }
   
  Double_t retval = sum/((Double_t)fN);
  if ( GetDebug() )
    {
      Info("AvarageDistance","Avarage distance is %f.",retval);
    }
  return retval;
}
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/

#include "AliRun.h"
#include "AliESD.h"
#include "AliRunLoader.h"
#include "AliTPCtrack.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "AliMagF.h"




void AliHBTTrackPoints::testesd(Int_t entr,const char* fname )
{
  delete gAlice;
  gAlice = 0x0;
  AliRunLoader* rl = AliRunLoader::Open();
  rl->LoadgAlice();
  
  Float_t mf = rl->GetAliRun()->Field()->SolenoidField();
 
  
  TFile* fFile = TFile::Open(fname);
  
  if (fFile == 0x0)
   {
     printf("testesd: There is no suche a ESD file\n");
     return;
   }
  AliESD* esd = dynamic_cast<AliESD*>(fFile->Get("0"));
  AliESDtrack *t = esd->GetTrack(entr);
  if (t == 0x0)
   {
     ::Error("testesd","Can not get track %d",entr);
     return;
   }

  
  Int_t N = 170;
  AliHBTTrackPoints* tp = new AliHBTTrackPoints(N,t,mf,1.);
  
  Float_t xmin = -250;
  Float_t xmax =  250;
  
  Float_t ymin = -250;
  Float_t ymax =  250;

  Float_t zmin = -250;
  Float_t zmax =  250;
  
  TH2D* hxy = new TH2D("hxy","hxy",1000,xmin,xmax,1000,ymin,ymax);
  TH2D* hxyt = new TH2D("hxyt","hxyt",1000,xmin,xmax,1000,ymin,ymax);
  TH2D* hxyTR = new TH2D("hxyTR","hxyTR",1000,xmin,xmax,1000,ymin,ymax);

  TH2D* hxz = new TH2D("hxz","hxz",1000,xmin,xmax,1000,zmin,zmax);
  TH2D* hxzt = new TH2D("hxzt","hxzt",1000,xmin,xmax,1000,zmin,zmax);
  TH2D* hxzTR = new TH2D("hxzTR","hxzTR",1000,xmin,xmax,1000,zmin,zmax);

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
      
   }

  rl->LoadTrackRefs();
  TTree* treeTR = rl->TreeTR();
  TBranch* b = treeTR->GetBranch("TPC");
  
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
  hxy->Draw("");
//  hxzt->Draw("same");
  hxyTR->Draw("same");
  
  delete rl;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/

void AliHBTTrackPoints::testtpc(Int_t entr)
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
  
  Float_t xmin = -250;
  Float_t xmax =  250;
  
  Float_t ymin = -250;
  Float_t ymax =  250;

  Float_t zmin = -250;
  Float_t zmax =  250;
  
  TH2D* hxy = new TH2D("hxy","hxy",1000,xmin,xmax,1000,ymin,ymax);
  TH2D* hxyt = new TH2D("hxyt","hxyt",1000,xmin,xmax,1000,ymin,ymax);
  TH2D* hxyTR = new TH2D("hxyTR","hxyTR",1000,xmin,xmax,1000,ymin,ymax);

  TH2D* hxz = new TH2D("hxz","hxz",1000,xmin,xmax,1000,zmin,zmax);
  TH2D* hxzt = new TH2D("hxzt","hxzt",1000,xmin,xmax,1000,zmin,zmax);
  TH2D* hxzTR = new TH2D("hxzTR","hxzTR",1000,xmin,xmax,1000,zmin,zmax);

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

