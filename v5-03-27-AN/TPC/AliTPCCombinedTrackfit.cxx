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
//
// Combine cosmic track pairs (upper, lower) and do track fitting
// oooooOOOOOooooo
// oooooOOOOOooooo
// oooooOOOOOooooo
//
//  Xianguo Lu <lu@physi.uni-heidelberg.de>

#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTreeStream.h>
#include <TVector3.h>

#include "AliESDtrack.h"
#include "AliESDfriendTrack.h"
#include "AliTPCseed.h"
#include "AliTrackerBase.h"
#include "AliTrackPointArray.h"

#include "AliTPCCosmicUtils.h"
#include "AliTPCCombinedTrackfit.h"

AliTPCCombinedTrackfit::AliTPCCombinedTrackfit(const Int_t dlev, const TString tag):
  fStreamer(0x0), fDebugLevel(dlev)
  , fSeedUp(0x0), fSeedLow(0x0), fTrackparUp(0x0), fTrackparLow(0x0)
  , fStatus(-999)
  , fLeverArm(-999)
  , fFitNcls(-999), fMissNcls(-999), fPreChi2(-999)
{
  //
  //Constructor
  //
  fInnerClusterUp.SetXYZ(-999,-999,-999);
  fInnerClusterLow.SetXYZ(-999,-999,-999);

  if(fDebugLevel>0)
    fStreamer = new TTreeSRedirector(Form("CombinedTrackfit_%s.root", tag.Data()));
}

AliTPCCombinedTrackfit::~AliTPCCombinedTrackfit()
{
  //
  //Destructor
  //

  delete fStreamer;
  
  delete fTrackparUp;
  delete fTrackparLow;
}

Bool_t AliTPCCombinedTrackfit::CombineESDtracks(AliESDtrack * &trk0, AliESDtrack *&trk1)
{
  //
  //Get TPCseeds from the 2 ESDtracks, swap TPCseeds and ESDTracks (if necessary) according to y (0:upper 1:lower), perform trackfit using TPCseeds
  //if fStatus==0, i.e. combine is successful, swap of the ESDtracks is kept since pointer *& is used
  //

  IniCombineESDtracks();

  if(!GetTPCseeds(trk0, trk1)){
    return kFALSE; 
  }

  Bool_t kswap = kFALSE;
  CombineTPCseeds(kswap);

  if(fStatus == 0){
    if(kswap){
      AliESDtrack * tmptrk = trk0;
      trk0 = trk1;
      trk1 = tmptrk;
    }
    return kTRUE;
  }
  else
    return kFALSE;
}

Bool_t AliTPCCombinedTrackfit::CombineTPCseeds(AliTPCseed * &seed0, AliTPCseed *&seed1)
{
  //
  //same as AliTPCCombinedTrackfit::CombineESDtracks, except that the seeds are passed in from outside, which can be still unordered
  //if fStatus==0, i.e. combine is successful, swap of the TPCseeds is kept since pointer *& is used
  //
  IniCombineESDtracks();

  fSeedUp  = seed0;
  fSeedLow = seed1;
  
  Bool_t kswap = kFALSE;
  CombineTPCseeds(kswap);

  if(fStatus==0){
    if(kswap){
      AliTPCseed * tmpseed = seed0;
      seed0 = seed1;
      seed1 = tmpseed;
    }
    return kTRUE;
  }
  else 
    return kFALSE;
}

void AliTPCCombinedTrackfit::Print() const
{
  //
  //print out variable values
  //
  printf("Status %2d NclsU %3d NclsD %3d ZinnerU %7.2f ZinnerD %7.2f LeverArm %7.2f\n", fStatus, fSeedUp->GetNumberOfClusters(), fSeedLow->GetNumberOfClusters(), fInnerClusterUp.Z(), fInnerClusterLow.Z(), fLeverArm);
}

Double_t AliTPCCombinedTrackfit::ImpactParameter() const
{
  //
  //calculate the impactparameter from (0,0,0)
  //
  const TVector3 p0(0,0,0);
  const TVector3 va = p0 - fInnerClusterUp;
  const TVector3 vb = fInnerClusterLow - fInnerClusterUp;

  const TVector3 dd = va.Cross(vb);

  return dd.Mag()/vb.Mag();
}

Double_t AliTPCCombinedTrackfit::MinPhi() const
{
  //
  //the smaller phi of the two tracks w.r.t. horizon
  //
  Double_t fsp[] = {fabs(sin(fTrackparUp->Phi())), fabs(sin(fTrackparLow->Phi()))};;
  return asin(TMath::Min(fsp[0], fsp[1])) * TMath::RadToDeg();
}
//===================================================================================================
//===================================================================================================

void AliTPCCombinedTrackfit::IniCombineESDtracks()
{
  //
  //initialization, for reuse of the same AliTPCCombinedTrackfit instance
  //

  fSeedUp = 0x0;
  fSeedLow = 0x0;
  delete fTrackparUp;
  delete fTrackparLow;
  fTrackparUp = 0x0;
  fTrackparLow = 0x0;

  fStatus = 0;
}

void AliTPCCombinedTrackfit::CombineTPCseeds(Bool_t &kswap)
{
  //
  //do combined trackfit using TPCseeds
  //

  if(
     !CheckNcls()
     || !AnaSeeds(kswap) 
     || !CheckLeverArm()
     )
    return;

  //AliExternalTrackParam object created
  fTrackparUp  = AliTPCCosmicUtils::MakeSeed(fSeedUp);
  fTrackparLow = AliTPCCosmicUtils::MakeSeed(fSeedLow);
  if(!fTrackparUp || !fTrackparLow){
    fStatus = kFailMakeSeed;
    return;
  }

  AliExternalTrackParam * trackPars[]={fTrackparUp, fTrackparLow};
  const AliTPCseed *seeds[]={fSeedUp, fSeedLow};
  TTreeSRedirector * debugstreamer = 0x0;
  if(fDebugLevel&2){
    debugstreamer = fStreamer;
  }

  AliTPCCosmicUtils::CombinedFit(trackPars, seeds, fFitNcls, fMissNcls, fPreChi2, debugstreamer);

  Update();

  return;
}

void AliTPCCombinedTrackfit::Update()
{
  //
  //Update variables depending on the fit result
  //

  if(fMissNcls || fFitNcls==0){
    fStatus = kFailPropagation;
    return;
  }

  fPreChi2 /= fFitNcls;
  if(fPreChi2>fgkMaxChi2){
    fStatus = kFailChi2;
    return;
  }

  if( fStatus == 0 && (fDebugLevel&1) ){
    Double_t momup  = fTrackparUp->P();
    Double_t momlow = fTrackparLow->P();
    Double_t ptup   = fTrackparUp->Pt();
    Double_t ptlow  = fTrackparLow->Pt();

    (*fStreamer)<<"TrackProp"<<
      "Tup.="<<fTrackparUp<<
      "Tlow.="<<fTrackparLow<<
      "icup.="<<&fInnerClusterUp<<
      "iclow.="<<&fInnerClusterLow<<
      "leverarm="<<fLeverArm<<
      "ncl="<<fFitNcls<<
      "nmiss="<<fMissNcls<<
      "chi2="<<fPreChi2<<
      "momup="<<  momup <<
      "momlow="<< momlow <<
      "ptup="<<   ptup <<
      "ptlow="<<  ptlow <<
      "\n";
  }
}

Bool_t AliTPCCombinedTrackfit::CheckLeverArm()
{
  //
  //if lever arm is too short, no need to use combined track fit. 
  //On the other hand, short lever arm from two tracks mostly means they are fake pairs.
  //lever arm extents over one quadrant, e.g. (0,250)-(250,0): 250*sqrt(2)~350
  //
  if(fLeverArm<fgkCutLeverArm){
    fStatus = kFailLeverArm;
    return kFALSE;
  }
  else 
    return kTRUE;
}

Bool_t AliTPCCombinedTrackfit::AnaSeeds(Bool_t &kswap)
{
  //
  //swap seeds (if necessary) so that (y of fSeedUp) > (y of fSeedLow)
  //

  //---------------------------------- navigate through all clusters ----------------------------------
  AliTPCseed ** seeds[]={&fSeedUp, &fSeedLow};

  //min, max according to y
  TVector3 singlemin[2], singlemax[2];
  for(Int_t ii=0; ii<2; ii++){
    singlemin[ii].SetXYZ( 1e10,  1e10,  1e10);
    singlemax[ii].SetXYZ(-1e10, -1e10, -1e10);
  }
  
  for(Int_t itrk=0; itrk<2; itrk++){
    for(Int_t irow=0; irow<AliTPCCosmicUtils::fgkNRow; irow++){
      const AliTPCclusterMI * cls = (*seeds[itrk])->GetClusterPointer(irow);
      if(!cls)
        continue;
      
      Float_t xyz[3]={-999,-999,-999};
      cls->GetGlobalXYZ(xyz);
      if(xyz[1]<singlemin[itrk].Y()){
        singlemin[itrk].SetXYZ(xyz[0], xyz[1], xyz[2]);
      }
      if(xyz[1]>singlemax[itrk].Y()){
        singlemax[itrk].SetXYZ(xyz[0], xyz[1], xyz[2]);
      }
    }
  }

  //--------------------------------

  //kpass true if y of the two seeds clearly separate: min of one > max of the other
  Bool_t kpass = kFALSE;

  fInnerClusterUp.SetXYZ(-999,-999,-999);
  fInnerClusterLow.SetXYZ(-999,-999,-999);
  TVector3 combinedmin, combinedmax;
  if(singlemin[0].Y()> singlemax[1].Y()){
    fInnerClusterUp  = singlemin[0];
    fInnerClusterLow = singlemax[1];

    //no need to swap
    kswap = kFALSE;

    kpass = kTRUE;

    combinedmax = singlemax[0];
    combinedmin = singlemin[1];
  }
  else if(singlemin[1].Y()> singlemax[0].Y()){
    fInnerClusterUp  = singlemin[1];
    fInnerClusterLow = singlemax[0];
  
    //have to be swapped
    kswap = kTRUE;
    AliTPCseed *tmp=*(seeds[0]);
    *(seeds[0])=*(seeds[1]);
    *(seeds[1])=tmp;
    
    kpass = kTRUE;

    combinedmax = singlemax[1];
    combinedmin = singlemin[0];
  }           
  else
    kpass = kFALSE;

  const TVector3 comdelta = combinedmax-combinedmin;
  fLeverArm = comdelta.Pt();

  if(!kpass){
    fStatus = kFailSwapSeeds;
    return kFALSE;
  }
  else
    return kTRUE;
}

Bool_t AliTPCCombinedTrackfit::CheckNcls()
{
  //
  //check number of clusters in TPCseed, for too small number MakeSeed will fail
  //
  if( fSeedUp->GetNumberOfClusters()<AliTPCCosmicUtils::fgkNclsMin || fSeedLow->GetNumberOfClusters()<AliTPCCosmicUtils::fgkNclsMin ){
    fStatus = kFailNclsMin;
    return kFALSE;
  }
  else
    return kTRUE;
}

Bool_t AliTPCCombinedTrackfit::GetTPCseeds(const AliESDtrack *trk0,  const AliESDtrack *trk1)
{
  //
  //Get TPC seeds from ESDfriendTrack
  //
  fSeedUp  = AliTPCCosmicUtils::GetTPCseed(trk0);
  fSeedLow = AliTPCCosmicUtils::GetTPCseed(trk1);

  if(!fSeedUp || !fSeedLow){
    fStatus = kFailGetTPCseeds;
    return kFALSE;
  }

  return kTRUE;
}

