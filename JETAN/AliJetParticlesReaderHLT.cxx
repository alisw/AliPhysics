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
#include <AliKalmanTrack.h>
#include <AliJetEventParticles.h>
#include "AliJetParticlesReaderHLT.h"

ClassImp(AliJetParticlesReaderHLT)

AliJetParticlesReaderHLT::AliJetParticlesReaderHLT(Bool_t bMapper, const Char_t* esdfilename) :
  AliJetParticlesReaderESD(esdfilename),
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
  AliJetParticlesReaderESD(dirs,esdfilename),
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

  Info("ReadESD","Reading Event %d",fCurrentEvent);
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
    if((fMinHits>0) && (knhits<fMinHits)) continue;
    if((fMinWeight>0) && (kweight<fMinWeight)) continue;

    const Float_t kpx=kesdtrack->GetPx();
    const Float_t kpy=kesdtrack->GetPy();
    const Float_t kpz=kesdtrack->GetPz();
    const Float_t kpt=kesdtrack->GetPt();
    const Float_t kp=TMath::Sqrt(kpz*kpz+kpt*kpt);
    const Float_t keta=0.5*TMath::Log((kp+kpz+1e-30)/(kp-kpz+1e-30)); 
    const Float_t kphi=TMath::Pi()+TMath::ATan2(-kpy,-kpx);

    if(IsAcceptedParticle(kpt,kphi,keta))
      fEventParticles->AddParticle(kpx,kpy,kpz,kp,i,kesdtrack->GetMCid(),kpt,kphi,keta);
  } //loop over tracks

  return kTRUE;
}
