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
#include <AliJetEventParticles.h>
#include "AliJetParticlesReaderHLT.h"

ClassImp(AliJetParticlesReaderHLT)

AliJetParticlesReaderHLT::AliJetParticlesReaderHLT(Bool_t bMapper, const Char_t* esdfilename) :
  AliJetParticlesReaderESD(esdfilename),
  fTrackerType(bMapper)
{
  //constructor
}

/********************************************************************/
  
AliJetParticlesReaderHLT::AliJetParticlesReaderHLT(
				      Bool_t bMapper,
                                      TObjArray* dirs,
                                      const Char_t* esdfilename) :
  AliJetParticlesReaderESD(dirs,esdfilename),
  fTrackerType(bMapper)
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

  Info("ReadESD","Reading Event %d",fCurrentEvent);
  if((!fOwner) || (fEventParticles==0)) 
    fEventParticles = new AliJetEventParticles();

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

  if(fTrackerType){
    const Int_t kntr =esd->GetNumberOfHLTHoughTracks();
    Info("ReadESD","Found %d conformal tracks.",kntr);
    for (Int_t i = 0;i<kntr; i++) {
      const AliESDHLTtrack *kesdtrack=esd->GetHLTConfMapTrack(i);

     if (kesdtrack == 0)
      {
        Error("ReadESD","Can not get track %d", i);
        continue;
      }

     //const Float_t kpid=kesdtrack->GetPID();
     //const Int_t knhits=kesdtrack->GetNHits();
     const Float_t kpx=kesdtrack->GetPx();
     const Float_t kpy=kesdtrack->GetPy();
     const Float_t kpz=kesdtrack->GetPz();
     const Float_t kpt=kesdtrack->GetPt();
     const Float_t kp=TMath::Sqrt(kpz*kpz+kpt*kpt);
     const Float_t keta=0.5*TMath::Log((kp+kpz+1e-30)/(kp-kpz+1e-30)); 
     const Float_t kphi=TMath::Pi()+TMath::ATan2(-kpy,-kpx);

     if(IsAcceptedParticle(kpt,kphi,keta))
       fEventParticles->AddParticle(kpx,kpy,kpz,kp,i,kesdtrack->GetMCid(),kpt,kphi,keta);

    } //loop over conf tracks
  } else {
    const Int_t kntr=esd->GetNumberOfHLTHoughTracks();
    Info("ReadESD","Found %d hough tracks.",kntr);

    for (Int_t i = 0;i<kntr; i++) {
      const AliESDHLTtrack *kesdtrack=esd->GetHLTHoughTrack(i);

     if (kesdtrack == 0)
      {
        Error("ReadESD","Can not get track %d", i);
        continue;
      }

     //const Float_t kweight=kesdtrack->GetWeight();
     const Float_t kpx=kesdtrack->GetPx();
     const Float_t kpy=kesdtrack->GetPy();
     const Float_t kpz=kesdtrack->GetPz();

     const Float_t kpt=kesdtrack->GetPt();
     const Float_t kp=TMath::Sqrt(kpz*kpz+kpt*kpt);
     const Float_t keta=0.5*TMath::Log((kp+kpz+1e-30)/(kp-kpz+1e-30)); 
     const Float_t kphi=TMath::Pi()+TMath::ATan2(-kpy,-kpx);

     if(IsAcceptedParticle(kpt,kphi,keta))
       fEventParticles->AddParticle(kpx,kpy,kpz,kp,i,kesdtrack->GetMCid(),kpt,kphi,keta);
    } //loop over hough tracks
  }
  return kTRUE;
}
