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


#include <TObjArray.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliTPCLaserTrack.h"

ClassImp(AliTPCLaserTrack)

TObjArray *AliTPCLaserTrack::fgArrLaserTracks=0x0;

AliTPCLaserTrack::AliTPCLaserTrack() :
  AliExternalTrackParam(),
  fId(-1),
  fSide(-1),
  fRod(-1),
  fBundle(-1),
  fBeam(-1),
  fRayLength(0)
{
  //
  // Default constructor
  //

}

AliTPCLaserTrack::AliTPCLaserTrack(const AliTPCLaserTrack &ltr) :
  AliExternalTrackParam(ltr),
  fId(ltr.fId),
  fSide(ltr.fSide),
  fRod(ltr.fRod),
  fBundle(ltr.fBundle),
  fBeam(ltr.fBeam),
  fRayLength(ltr.fRayLength)
{
  //
  // Default constructor
  //

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
  fRayLength(rayLength)
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

void AliTPCLaserTrack::LoadTracks()
{
    //
    // Load all design positions from file into the static array fgArrLaserTracks
    //

    if ( fgArrLaserTracks ) return;

    TString dataFileName("$ALICE_ROOT/TPC/Calib/LaserTracks.root");  //Path to the Data File

    TFile *f=TFile::Open(gSystem->ExpandPathName(dataFileName.Data()));
    if ( !f || !f->IsOpen() ){
//	AliWarning(Form("Could not open laser data file: '%s'",dataFileName.Data()));
//	AliWarning("Could not open laser data file");
	return;
    }
    TObjArray *arrLaserTracks = (TObjArray*)f->Get("arrLaserTracks");
    if ( !arrLaserTracks ) {
//	AliWarning(Form("Could not get laser position data from file: '%s'",fgkDataFileName));
        return;
    }

    fgArrLaserTracks = new TObjArray(fgkNLaserTracks);
    for (Int_t itrack=0; itrack<fgkNLaserTracks; itrack++){
	AliTPCLaserTrack *ltr = (AliTPCLaserTrack*)arrLaserTracks->At(itrack);
	if ( !ltr ){
//	    AliWarning(Form("No informatino found for Track %d!",itrack));
	    continue;
	}
        fgArrLaserTracks->AddAt(new AliTPCLaserTrack(*ltr),itrack);
    }
    delete f;
}


Int_t AliTPCLaserTrack::IdentifyTrack(AliExternalTrackParam *track)
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
  track->GetPxPyPz(pxyz0);
  //
  Float_t mindist=40; // maxima minimal distance
  Int_t id = -1;
  AliExternalTrackParam*  ltr0= (AliExternalTrackParam*)arrTracks->UncheckedAt(0);
  for (Int_t itrack=0; itrack<fgkNLaserTracks; itrack++){
    AliExternalTrackParam *ltr = (AliExternalTrackParam*)arrTracks->UncheckedAt(itrack);
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
    ltr->GetPxPyPz(pxyz1);
    Float_t distdir = (ltr->GetParameter()[2]-track->GetParameter()[2])*90; //distance at entrance
    if (TMath::Abs(ltr->GetParameter()[2]-track->GetParameter()[2])>kMaxdphiP)
      continue;
    //
    Float_t dist=0;
    dist+=TMath::Abs(lxyz1[0]-lxyz0[0]);
    dist+=TMath::Abs(lxyz1[1]-lxyz0[1]);
    dist+=TMath::Abs(lxyz1[2]-lxyz0[2]);
    dist+=distdir;
    //    
    if (id<0)  {
      id =itrack; 
      mindist=dist; 
      ltr0=ltr;
      continue;
    }
    if (dist>mindist) continue;
    id = itrack;
    mindist=dist;
    ltr0=ltr;
  }
  return id;
}

