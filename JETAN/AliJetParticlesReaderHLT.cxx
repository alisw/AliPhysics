// $Id$

//___________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// File reader for HLT tracks ESD                                          //
//                                                                         //
// loizides@ikf.uni-frankfurt.de                                           //
/////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>
#include <TString.h>
#include <TObjString.h>
#include <TTree.h>
#include <TFile.h>
#include <AliESD.h>
#include <AliESDtrack.h>
#include <AliESDHLTtrack.h>
#include <AliL3Track.h>
#include <AliL3Vertex.h>
#include <AliKalmanTrack.h>
#include <AliJetEventParticles.h>
#include "AliJetParticlesReaderHLT.h"

ClassImp(AliJetParticlesReaderHLT)

AliJetParticlesReaderHLT::AliJetParticlesReaderHLT(Bool_t bMapper, const Char_t* esdfilename) :
  AliJetParticlesReaderESD(0,esdfilename),
  fTrackerType(bMapper),
  fMinHits(0),
  fMinWeight(0)
{
  //constructor
}

/********************************************************************/
  
AliJetParticlesReaderHLT::AliJetParticlesReaderHLT(
				      Bool_t bMapper,
                                      TObjArray* dirs,
                                      const Char_t* esdfilename) :
  AliJetParticlesReaderESD(0,dirs,esdfilename),
  fTrackerType(bMapper),
  fMinHits(0),
  fMinWeight(0)

{
  //constructor
}

/********************************************************************/

AliJetParticlesReaderHLT::~AliJetParticlesReaderHLT()
{
  //desctructor
}

Int_t AliJetParticlesReaderHLT::ReadESD(AliESD* esd)
{
  //Reads one ESD

  if (esd == 0)
   {
     Error("ReadESD","ESD is NULL");
     return kFALSE;
   }

  Float_t mf = esd->GetMagneticField(); 
  if (mf <= 0.0)  
  {
     Error("ReadESD","Magnetic Field is 0. Skipping to next event.");
     return kFALSE;
  }
  AliKalmanTrack::SetMagneticField(mf/10.);

  Info("ReadESD","Reading Event %d",fCurrentDir*1000+fCurrentEvent);
  if((!fOwner) || (fEventParticles==0)) 
    fEventParticles = new AliJetEventParticles();

  Int_t ntr=0;
  if(fTrackerType){
    ntr =esd->GetNumberOfHLTHoughTracks();
    Info("ReadESD","Found %d conformal tracks.",ntr);
  } else {
    ntr=esd->GetNumberOfHLTHoughTracks();
    Info("ReadESD","Found %d hough tracks.",ntr);
  }
  fEventParticles->Reset(ntr);

  TString headdesc="";
  headdesc+="Run ";
  headdesc+=esd->GetRunNumber();
  headdesc+=": Ev ";
  headdesc+=esd->GetEventNumber();
  fEventParticles->SetHeader(headdesc);

  Double_t vertexpos[3];//vertex position
  const AliESDVertex* kvertex = esd->GetVertex();
  if (kvertex == 0)
   {
     Info("ReadESD","ESD returned NULL pointer to vertex - assuming (0.0,0.0,0.0)");
     vertexpos[0] = 0.0;
     vertexpos[1] = 0.0;
     vertexpos[2] = 0.0;
   }
  else
   {
     kvertex->GetXYZ(vertexpos);
   }

  fEventParticles->SetVertex(vertexpos[0],vertexpos[1],vertexpos[2]);
  //cout << vertexpos[0] << " " << vertexpos[1] << " " << vertexpos[2] << endl;

  AliL3Track l3;
  AliL3Vertex v;
  v.SetX(vertexpos[0]);
  v.SetY(vertexpos[1]);
  v.SetZ(vertexpos[2]);
  
  Double_t xc=0.,yc=0.,zc=0.;
  for (Int_t i = 0;i<ntr; i++) {
    AliESDHLTtrack *kesdtrack;
    if(fTrackerType){
      kesdtrack=esd->GetHLTConfMapTrack(i);
    } else {
      kesdtrack=esd->GetHLTHoughTrack(i);
    }

    if (kesdtrack == 0)
      {
        Error("ReadESD","Can not get track %d", i);
        continue;
      }

    //const Float_t kpid=kesdtrack->GetPID();
    const Int_t knhits=kesdtrack->GetNHits();
    const Int_t kweight=kesdtrack->GetWeight();
    //cout << i << " " << kweight << " " << knhits << endl;
    if((fMinHits>0) && (knhits<fMinHits)) continue;    
    if(kweight>1000) continue; //avoid ghosts 
    if((fMinWeight>0) && (kweight<fMinWeight)) continue;

    Float_t px=kesdtrack->GetPx();
    Float_t py=kesdtrack->GetPy();
    Float_t pz=kesdtrack->GetPz();

    if(0&&fTrackerType){
      //if(!kesdtrack->ComesFromMainVertex()) continue;
      cout << kesdtrack->GetPx() << " " << kesdtrack->GetPy() << " " << kesdtrack->GetPt() << endl;
      l3.SetFirstPoint(kesdtrack->GetFirstPointX(),kesdtrack->GetFirstPointY(),kesdtrack->GetFirstPointZ());
      l3.SetLastPoint(kesdtrack->GetLastPointX(),kesdtrack->GetLastPointY(),kesdtrack->GetLastPointZ());
      l3.SetCharge(kesdtrack->GetCharge());
      l3.SetPt(kesdtrack->GetPt());
      l3.SetTgl(kesdtrack->GetTgl());
      l3.SetPsi(kesdtrack->GetPsi());
      l3.CalculateHelix();
      l3.GetClosestPoint(&v,xc,yc,zc);
      if(TMath::Abs(zc)>10.) continue;
      l3.SetFirstPoint(vertexpos[0],vertexpos[1],vertexpos[2]);
      //l3.CalculateHelix();
      l3.UpdateToFirstPoint();
      px=l3.GetPx();
      py=l3.GetPy();
      pz=l3.GetPz();
    }

    const Float_t kpt=TMath::Sqrt(px*px+py*py);
    //cout << px << " " << py << " " << " " << kpt << endl;

    const Float_t kp=TMath::Sqrt(pz*pz+kpt*kpt);
    const Float_t keta=0.5*TMath::Log((kp+pz+1e-30)/(kp-pz+1e-30)); 
    const Float_t kphi=TMath::Pi()+TMath::ATan2(-py,-px);

    if(IsAcceptedParticle(kpt,kphi,keta))
      fEventParticles->AddParticle(px,py,pz,kp,i,kesdtrack->GetMCid(),knhits,kpt,kphi,keta);
  } //loop over tracks

  return kTRUE;
}

#if 0
  SetChi2(0.);
  if(t.GetNHits()==1)
    SetNumberOfClusters(0);
  else
    SetNumberOfClusters(t.GetNHits());
  SetLabel(t.GetMCid());
  SetMass(0.13957);

  fdEdx=0;
  fAlpha = fmod((t.GetSector()+0.5)*(2*TMath::Pi()/18),2*TMath::Pi());
  if      (fAlpha < -TMath::Pi()) fAlpha += 2*TMath::Pi();
  else if (fAlpha >= TMath::Pi()) fAlpha -= 2*TMath::Pi();

  //First the emiision angle
  Double_t psi = t.GetPsi()-(t.GetSector()+0.5)*(2*TMath::Pi()/18);

  //Then local x,y coordinates
  Double_t radius = t.GetPt()*GetConvConst();
  Double_t xhit = 82.97; //Position at first TPC padrow
  Double_t trackphi0 = psi + (-t.GetCharge())*TMath::Pi()/2;
  Double_t x0 = t.GetFirstPointX()*TMath::Cos(fAlpha) + t.GetFirstPointY()*TMath::Sin(fAlpha);
  Double_t y0 = t.GetFirstPointY()*TMath::Cos(fAlpha) - t.GetFirstPointX()*TMath::Sin(fAlpha);
  Double_t centerx = radius *  cos(trackphi0) + x0;
  Double_t centery = radius *  sin(trackphi0) + y0;
  Double_t aa = (xhit - centerx)*(xhit - centerx);
  Double_t r2 = radius*radius;
  if(aa > r2) throw "AliITStrackV2: conversion failed !\n";
  Double_t aa2 = sqrt(r2 - aa);
  Double_t y1 = centery + aa2;
  Double_t y2 = centery - aa2;
  Double_t yhit = y1;
  if(fabs(y2) < fabs(y1)) yhit = y2;

  //Local z coordinate
  Double_t angle1 = atan2((yhit - centery),(xhit - centerx));
  if(angle1 < 0) angle1 += 2.*TMath::Pi();
  Double_t angle2 = atan2((x0-centery),(y0-centerx));
  if(angle2 < 0) angle2 += 2.*TMath::Pi();
  Double_t diffangle = angle1 - angle2;
  diffangle = fmod(diffangle,2.*TMath::Pi());
  if(((-t.GetCharge())*diffangle) < 0) diffangle = diffangle - (-t.GetCharge())*2.*TMath::Pi();
  Double_t stot = fabs(diffangle)*radius;
  Double_t zhit;
  if(t.GetNHits()==1)
    zhit = zvertex + stot*t.GetTgl();
  else
    zhit = t.GetFirstPointZ() + stot*t.GetTgl();

  //Local sine of track azimuthal angle
  if((-t.GetCharge())<0)
    radius = -radius;
  Double_t sinbeta = -1.*(centerx - xhit)/radius;

  //Filling of the track paramaters
  fX=xhit;
  fP0=yhit;
  fP1=zhit;
  fP2=sinbeta;
  fP3=t.GetTgl();
  fP4=1./radius;
#endif

