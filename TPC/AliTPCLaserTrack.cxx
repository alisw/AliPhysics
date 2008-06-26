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

    LoadTracks();
    TObjArray *arrTracks = GetTracks();

    Double_t phitr=0;
    Double_t philtr=0;
    Double_t xltr[3];
    Double_t xtr[3];
    Double_t vtr[3];

    track->GetXYZ(xtr);
    track->GetDirection(vtr);
    phitr=track->Phi();

    Int_t id = -1;
    Double_t dcaMin=3;  // 3 sigma is the minimum weighted dca accepted
    for (Int_t itrack=0; itrack<fgkNLaserTracks; itrack++){
	AliExternalTrackParam *ltr = (AliExternalTrackParam*)arrTracks->UncheckedAt(itrack);
	philtr=ltr->Phi();
	Double_t phiadd=0;
        Double_t dphi=TMath::Abs(phitr-philtr);
	if (dphi>TMath::Pi()) phiadd=TMath::Pi();
        dphi-=phiadd;
	if (dphi>2*TMath::DegToRad()) continue;  //only 2 degree in phi

	//printf("itrack: %d; dphi: %f\n",itrack,dphi/TMath::DegToRad());

	ltr->GetXYZ(xltr);

	Double_t l=0;
	Double_t d2=0;
	Double_t d=0;
        for (Int_t i=0; i<2; i++)
	    l+=(xltr[i]-xtr[i])*vtr[i];
        for (Int_t i=0; i<2; i++)
	    d2+=(xltr[i]-xtr[i]-l*vtr[i])*(xltr[i]-xtr[i]-l*vtr[i]);
	d=TMath::Sqrt(d2);

	//printf("itrack: %d; d-xy: %f\n",itrack,d);


        //only 3mm in x-y
        if ( d>.3 ) continue;


	Double_t dz=TMath::Abs(xltr[2] - xtr[2]);
	//printf("itrack: %d; d-z: %f\n",itrack,dz);
        //30 cm in z
	if ( dz > 30. ) continue;

	//if ( id!=-1 ) printf("Warnig: Track (%d) already identified before (%d)\n", itrack, id);
        id=itrack;
//        Double_t relDistX
//        Double_t xtr=0;
//        Double_t xltr=0;
//	Double_t dca=track->GetDCA(ltr,0,xtr,xltr);
//	if ( dca<dcaMin ){
//	    id=itrack;
//	    dcaMin=dca;
//        }
    }
    return id;
}

