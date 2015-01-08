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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Surveyed Laser Track positions                                         //
// the position and direction information are stored in                   //
// the AliExternalTrackParam base class                                   //
// This class extends this information by identification parameters       //
/*

//Dump positions to a tree:
AliTPCLaserTrack::LoadTracks();
TObjArray *arr=AliTPCLaserTrack::GetTracks();
TTreeSRedirector *s=new TTreeSRedirector("LaserTracks.root");
TIter next(arr);
TObject *o=0x0;
while ( (o=next()) ) (*s) << "tracks" << "l.=" << o << "\n";
delete s;

//draw something
TFile f("LaserTracks.root");
TTree *tracks=(TTree*)f.Get("tracks");
tracks->Draw("fVecGY.fElements:fVecGX.fElements");

 tracks->Draw("fVecGY.fElements:fVecGX.fElements>>h(500,-250,250,500,-250,250)","fId<7")
*/
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <TObjArray.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#include "AliTPCLaserTrack.h"
#include "AliTPCROC.h"

ClassImp(AliTPCLaserTrack)

TObjArray *AliTPCLaserTrack::fgArrLaserTracks=0x0;

AliTPCLaserTrack::AliTPCLaserTrack() :
  AliExternalTrackParam(),
  fId(-1),
  fSide(-1),
  fRod(-1),
  fBundle(-1),
  fBeam(-1),
  fRayLength(0),
  fVecSec(0),       // points vectors - sector
  fVecP2(0),       // points vectors - snp
  fVecPhi(0),       // points vectors - global phi
  fVecGX(0),       // points vectors - globalX
  fVecGY(0),       // points vectors - globalY
  fVecGZ(0),       // points vectors - globalZ
  fVecLX(0),       // points vectors - localX
  fVecLY(0),       // points vectors - localY
  fVecLZ(0)        // points vectors - localZ
{
  //
//   // Default constructor
  //

}

AliTPCLaserTrack::AliTPCLaserTrack(const AliTPCLaserTrack &ltr) :
  AliExternalTrackParam(ltr),
  fId(ltr.fId),
  fSide(ltr.fSide),
  fRod(ltr.fRod),
  fBundle(ltr.fBundle),
  fBeam(ltr.fBeam),
  fRayLength(ltr.fRayLength),
  fVecSec(0),       // points vectors - sector
  fVecP2(0),       // points vectors - snp
  fVecPhi(0),       // points vectors - global phi
  fVecGX(0),       // points vectors - globalX
  fVecGY(0),       // points vectors - globalY
  fVecGZ(0),       // points vectors - globalZ
  fVecLX(0),       // points vectors - localX
  fVecLY(0),       // points vectors - localY
  fVecLZ(0)        // points vectors - localZ
{
  //
  // Default constructor
  //
  fVecSec=new TVectorD(*ltr.fVecSec);       // points vectors - sector
  fVecP2 =new TVectorD(*ltr.fVecP2);       // points vectors - snp
  fVecPhi=new TVectorD(*ltr.fVecPhi);       // points vectors - global phi
  fVecGX =new TVectorD(*ltr.fVecGX);       // points vectors - globalX
  fVecGY =new TVectorD(*ltr.fVecGY);       // points vectors - globalY
  fVecGZ =new TVectorD(*ltr.fVecGZ);       // points vectors - globalZ
  fVecLX =new TVectorD(*ltr.fVecLX);       // points vectors - localX
  fVecLY =new TVectorD(*ltr.fVecLY);       // points vectors - localY
  fVecLZ =new TVectorD(*ltr.fVecLZ);       // points vectors - localY

}

AliTPCLaserTrack::AliTPCLaserTrack(const Int_t id, const Int_t side, const Int_t rod,
		     const Int_t bundle, const Int_t beam,
		     Double_t x, Double_t alpha,
		     const Double_t param[5],
		     const Double_t covar[15], const Float_t rayLength) :
  AliExternalTrackParam(x,alpha,param,covar),
  fId(id),
  fSide(side),
  fRod(rod),
  fBundle(bundle),
  fBeam(beam),
  fRayLength(rayLength),
  fVecSec(new TVectorD(159)),       // points vectors - sector
  fVecP2(new TVectorD(159)),       // points vectors - snp
  fVecPhi(new TVectorD(159)),       // points vectors - global phi
  fVecGX(new TVectorD(159)),       // points vectors - globalX
  fVecGY(new TVectorD(159)),       // points vectors - globalY
  fVecGZ(new TVectorD(159)),       // points vectors - globalZ
  fVecLX(new TVectorD(159)),       // points vectors - localX
  fVecLY(new TVectorD(159)),       // points vectors - localY
  fVecLZ(new TVectorD(159))        // points vectors - localZ

{
  //
  // create laser track from arguments
  //
  
}
//_____________________________________________________________________
AliTPCLaserTrack& AliTPCLaserTrack::operator = (const  AliTPCLaserTrack &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTPCLaserTrack(source);
  
  return *this;
}


AliTPCLaserTrack::~AliTPCLaserTrack(){
  //
  // destructor
  //
  delete fVecSec;      //                - sector numbers  
  delete fVecP2;       //                - P2  
  delete fVecPhi;       // points vectors - global phi
  delete fVecGX;       // points vectors - globalX
  delete fVecGY;       // points vectors - globalY
  delete fVecGZ;       // points vectors - globalZ
  delete fVecLX;       // points vectors - localX
  delete fVecLY;       // points vectors - localY
  delete fVecLZ;       // points vectors - localZ
}

void AliTPCLaserTrack::LoadTracks()
{
  //
  // Load all design positions from file into the static array fgArrLaserTracks
  //
  
  if ( fgArrLaserTracks ) return;
  TObjArray *arrLaserTracks = 0x0;
  
  AliCDBManager *man=AliCDBManager::Instance();
  if (!man->GetDefaultStorage() && gSystem->Getenv("ALICE_ROOT")) man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  if (man->GetDefaultStorage()){
    if (man->GetRun()<0) man->SetRun(0);
    AliCDBEntry *entry=man->Get(AliCDBPath("TPC/Calib/LaserTracks"));
    if (!entry) return;
    arrLaserTracks = (TObjArray*)entry->GetObject();
    entry->SetOwner(kTRUE);
  } else {
    if (!gSystem->AccessPathName("LaserTracks.root")){
      TFile f("LaserTracks.root");
      arrLaserTracks=(TObjArray*)f.Get("arrLaserTracks");
      f.Close();
    }
  }
  if ( !arrLaserTracks ) {
//	AliWarning(Form("Could not get laser position data from file: '%s'",fgkDataFileName));
    return;
  }
  
  arrLaserTracks->SetOwner();
  
  fgArrLaserTracks = new TObjArray(fgkNLaserTracks);
  fgArrLaserTracks->SetOwner();
  for (Int_t itrack=0; itrack<fgkNLaserTracks; itrack++){
    AliTPCLaserTrack *ltr = (AliTPCLaserTrack*)arrLaserTracks->At(itrack);
    if ( !ltr ){
//	    AliWarning(Form("No informatino found for Track %d!",itrack));
      continue;
    }
    ltr->UpdatePoints();
    fgArrLaserTracks->AddAt(new AliTPCLaserTrack(*ltr),itrack);
  }

  delete arrLaserTracks;
}


void AliTPCLaserTrack::UpdatePoints(){
  //
  // update track points
  //
  const Double_t kMaxSnp=0.97;
  AliTPCROC* roc = AliTPCROC::Instance();
  //
  //
  if (!fVecSec){
    fVecSec=new TVectorD(159);
    fVecP2 =new TVectorD(159);       //                - P2  
    fVecPhi=new TVectorD(159);       //                - Phi
    fVecGX=new TVectorD(159);       // points vectors - globalX
    fVecGY=new TVectorD(159);       // points vectors - globalY
    fVecGZ=new TVectorD(159);       // points vectors - globalZ
    fVecLX=new TVectorD(159);       // points vectors - localX
    fVecLY=new TVectorD(159);       // points vectors - localY
    fVecLZ=new TVectorD(159);       // points vectors - localZ

  }
  for (Int_t irow=158; irow>=0; irow--){
    (*fVecSec)[irow]= -1;       //                -
    (*fVecP2)[irow] = 0;       //                - P2  -snp
    (*fVecPhi)[irow]= 0;       //                - global phi
    (*fVecGX)[irow] = 0;       // points vectors - globalX
    (*fVecGY)[irow] = 0;       // points vectors - globalY
    (*fVecGZ)[irow] = 0;       // points vectors - globalZ
    (*fVecLX)[irow] = 0;       // points vectors - localX
    (*fVecLY)[irow] = 0;       // points vectors - localY
    (*fVecLZ)[irow] = 0;       // points vectors - localZ

  }
  Double_t gxyz[3];
  Double_t lxyz[3];
  AliTPCLaserTrack*ltrp=new AliTPCLaserTrack(*this);  //make temporary track

  for (Int_t irow=158; irow>=0; irow--){
    UInt_t srow = irow;
    Int_t sector=0;
   
    if (srow >=roc->GetNRows(0)) {
      srow-=roc->GetNRows(0);
      sector=36    ;
    }
    lxyz[0]= roc->GetPadRowRadii(sector,srow);
    if (!ltrp->PropagateTo(lxyz[0],5)) break;
    ltrp->GetXYZ(gxyz);
    //
    Double_t alpha=TMath::ATan2(gxyz[1],gxyz[0]);
    if (alpha<0) alpha+=2*TMath::Pi();
    sector      +=TMath::Nint(-0.5+9*alpha/TMath::Pi());
    if (gxyz[2]<0) sector+=18;
    Double_t salpha   = TMath::Pi()*(sector+0.5)/9.;    
    if (!ltrp->Rotate(salpha)) break;
    if (!ltrp->PropagateTo(lxyz[0],5)) break;
    if (TMath::Abs(ltrp->GetSnp())>kMaxSnp) break;
    ltrp->GetXYZ(gxyz);
    lxyz[1]=ltrp->GetY();
    lxyz[2]=ltrp->GetZ();
    (*fVecSec)[irow]= sector;
    (*fVecP2)[irow] = ltrp->GetSnp();                 //                - P2  -snp
    (*fVecPhi)[irow]= TMath::ATan2(gxyz[1],gxyz[0]);  //                - global phi
    (*fVecGX)[irow] = gxyz[0];       // points vectors - globalX
    (*fVecGY)[irow] = gxyz[1];       // points vectors - globalY
    (*fVecGZ)[irow] = gxyz[2];       // points vectors - globalZ
    (*fVecLX)[irow] = lxyz[0];       // points vectors - localX
    (*fVecLY)[irow] = lxyz[1];       // points vectors - localY
    (*fVecLZ)[irow] = lxyz[2];       // points vectors - localZ

  }
  delete ltrp;  // delete temporary track
}

Int_t AliTPCLaserTrack::IdentifyTrack(AliExternalTrackParam *track, Int_t side)
{
  //
  // Find the laser track which is corresponding closest to 'track'
  // return its id
  //
  // 
  const  Float_t   kMaxdphi=0.2;
  const  Float_t   kMaxdphiP=0.05;
  const  Float_t   kMaxdz=40;

  if ( !fgArrLaserTracks ) LoadTracks();
  TObjArray *arrTracks = GetTracks();
  Double_t lxyz0[3];
  Double_t lxyz1[3];
  Double_t pxyz0[3];
  Double_t pxyz1[3];
  track->GetXYZ(lxyz0);
  track->GetDirection(pxyz0);
  //
  Float_t mindist=10; // maxima minimal distance
  Int_t id = -1;
  for (Int_t itrack=0; itrack<fgkNLaserTracks; itrack++){    
    AliTPCLaserTrack *ltr = (AliTPCLaserTrack*)arrTracks->UncheckedAt(itrack);
    if (side>=0) if (ltr->GetSide()!=side) continue;
    Double_t * kokot = (Double_t*)ltr->GetParameter();
    kokot[4]=-0.0000000001;
    //
    ltr->GetXYZ(lxyz1);
    if (TMath::Abs(lxyz1[2]-lxyz0[2])>kMaxdz) continue;
    // phi position
    Double_t phi0 = TMath::ATan2(lxyz0[1],lxyz0[0]);
    Double_t phi1 = TMath::ATan2(lxyz1[1],lxyz1[0]);
    if (TMath::Abs(phi0-phi1)>kMaxdphi) continue;
    // phi direction
    ltr->GetDirection(pxyz1);
    Float_t direction= pxyz0[0]*pxyz1[0] + pxyz0[1]*pxyz1[1] + pxyz0[2]*pxyz1[2];
    Float_t distdir = (1-TMath::Abs(direction))*90.; //distance at entrance
    if (1-TMath::Abs(direction)>kMaxdphiP)
      continue;
    //
    Float_t dist=0;
    dist+=TMath::Abs(lxyz1[0]-lxyz0[0]);
    dist+=TMath::Abs(lxyz1[1]-lxyz0[1]);
    //    dist+=TMath::Abs(lxyz1[2]-lxyz0[2]); //z is not used for distance calculation
    dist+=distdir;
    //    
    if (id<0)  {
      id =itrack; 
      mindist=dist; 
      continue;
    }
    if (dist>mindist) continue;
    id = itrack;
    mindist=dist;
  }
  return id;
}

