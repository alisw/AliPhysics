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
#include <AliHLTTrack.h>
#include <AliHLTVertex.h>
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
  AliHLTVertex v;
  v.SetX(vertexpos[0]);
  v.SetY(vertexpos[1]);
  v.SetZ(vertexpos[2]);

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
    if(TMath::IsNaN(px)||TMath::IsNaN(py)||TMath::IsNaN(pz)) continue;
    if(TMath::Abs(px)>1e3||TMath::Abs(py)>1e3||TMath::Abs(pz)>1e3) continue;
    //cout << px << " " << py << " " << pz <<  " " << TMath::Sqrt(px*px+py*py) << endl;

    if(0&&fTrackerType){
#if 0
      //kesdtrack->SetCharge(-kesdtrack->GetCharge());
	Double_t mom[3];
	if(!kesdtrack->GetPxPyPzAt(0,mom)) continue;
	px=mom[0];
	py=mom[1];
	pz=mom[2];
#else
      AliHLTTrack l3;
      //if(!kesdtrack->ComesFromMainVertex()) continue;
      //cout << "Pos: " << kesdtrack->GetFirstPointX() << " " << kesdtrack->GetFirstPointY() << " " << kesdtrack->GetFirstPointZ() << endl;      
      l3.SetFirstPoint(kesdtrack->GetFirstPointX(),kesdtrack->GetFirstPointY(),kesdtrack->GetFirstPointZ());
      l3.SetLastPoint(kesdtrack->GetLastPointX(),kesdtrack->GetLastPointY(),kesdtrack->GetLastPointZ());
      l3.SetCharge(kesdtrack->GetCharge());
      l3.SetPt(kesdtrack->GetPt());
      l3.SetTgl(kesdtrack->GetTgl());
      l3.SetPsi(kesdtrack->GetPsi());
      l3.CalculateHelix();
      Double_t xc=0.,yc=0.,zc=0.;
      l3.GetClosestPoint(&v,xc,yc,zc);
      if(TMath::Abs(zc)>10.) continue;
      l3.SetFirstPoint(0,0,0);
      //cout << "Pos: " << xc << " " << yc << " " << zc << endl;
      l3.UpdateToFirstPoint();
      px=l3.GetPx();
      py=l3.GetPy();
      pz=l3.GetPz();
#endif
    }
    const Float_t kpt=TMath::Sqrt(px*px+py*py);
    //cout << px << " " << py << " " << pz <<  " " << kpt << endl;

    const Float_t kp=TMath::Sqrt(pz*pz+kpt*kpt);
    const Float_t keta=0.5*TMath::Log((kp+pz+1e-30)/(kp-pz+1e-30)); 
    const Float_t kphi=TMath::Pi()+TMath::ATan2(-py,-px);

    if(IsAcceptedParticle(kpt,kphi,keta))
      fEventParticles->AddParticle(px,py,pz,kp,i,kesdtrack->GetMCid(),knhits,kpt,kphi,keta);
  } //loop over tracks

  return kTRUE;
}
