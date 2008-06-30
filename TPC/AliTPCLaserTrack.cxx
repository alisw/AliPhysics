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
  fBeam(-1)
{
  //
  // Default constructor
  //

}

AliTPCLaserTrack::AliTPCLaserTrack(AliTPCLaserTrack &ltr) :
  AliExternalTrackParam(ltr),
  fId(ltr.fId),
  fSide(ltr.fSide),
  fRod(ltr.fRod),
  fBundle(ltr.fBundle),
  fBeam(ltr.fBeam)
{
  //
  // Default constructor
  //

}

AliTPCLaserTrack::AliTPCLaserTrack(const Int_t id, const Int_t side, const Int_t rod,
		     const Int_t bundle, const Int_t beam,
		     Double_t x, Double_t alpha,
		     const Double_t param[5],
		     const Double_t covar[15]) :
  AliExternalTrackParam(x,alpha,param,covar),
  fId(id),
  fSide(side),
  fRod(rod),
  fBundle(bundle),
  fBeam(beam)
{
  //
  // create laser track from arguments
  //

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
  const  Float_t   kMaxdphi=0.1;
  const  Float_t   kMaxdphiP=0.06;
  const  Float_t   kMaxdz=50;

  if ( !fgArrLaserTracks ) LoadTracks();
  TObjArray *arrTracks = GetTracks();
  

  Double_t lxyz0[3];
  Double_t lxyz1[3];
  Double_t pxyz0[3];
  Double_t pxyz1[3];
  track->GetXYZ(lxyz0);
  track->GetPxPyPz(pxyz0);
  
  Int_t id = -1;
  for (Int_t itrack=0; itrack<fgkNLaserTracks; itrack++){
    AliExternalTrackParam *ltr = (AliExternalTrackParam*)arrTracks->UncheckedAt(itrack);
    Double_t * kokot = (Double_t*)ltr->GetParameter();
    kokot[4]=-0.0000000001;
    //
    ltr->GetXYZ(lxyz1);
    if ( (lxyz1[2]>0) && lxyz0[2]<0) continue;
    if ( (lxyz1[2]<0) && lxyz0[2]>0) continue;
    if (TMath::Abs(lxyz1[2]-lxyz0[2])>kMaxdz) continue;
    // phi position
    Double_t phi0 = TMath::ATan2(lxyz0[1],lxyz0[0]);
    Double_t phi1 = TMath::ATan2(lxyz1[1],lxyz1[0]);
    if (TMath::Abs(phi0-phi1)>kMaxdphi) continue;
    // phi direction
    ltr->GetPxPyPz(pxyz1);
    Double_t pphi0 = TMath::ATan2(pxyz0[1],pxyz0[0]);
    Double_t pphi1 = TMath::ATan2(pxyz1[1],pxyz1[0]);
    Bool_t phimatch = kFALSE;
    if (TMath::Abs(ltr->GetParameter()[2]-track->GetParameter()[2])>kMaxdphiP)
      continue;
  //   if (TMath::Abs(pphi0-pphi1)<kMaxdphiP) phimatch=kTRUE;
//     if (TMath::Abs(pphi0-pphi1-TMath::Pi())<kMaxdphiP) phimatch=kTRUE;
//     if (TMath::Abs(pphi0-pphi1+TMath::Pi())<kMaxdphiP) phimatch=kTRUE;
//     if (!phimatch) continue;
    //
    id =itrack;
  }
  return id;
}

