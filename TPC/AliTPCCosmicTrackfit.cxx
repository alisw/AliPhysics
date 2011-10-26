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
//-------------------------------------------------------
//          Implementation of the TPC combined trackfit
//
//   Origin:  Xianguo Lu lu@physi.uni-heidelberg.de  Xianguo.Lu@cern.ch
//
//=========================================================================================
// Motivation:
// 
// In the default reconstruction in the ALICE the cosmic tracks are found  as two independent particles.
//
// In general any of  subtracks can be used for the physics studies. In order to avoid the double counting,
// the track from the upper hemisphere can be used. 
//
// The momentum resolution is determined by the lever arm (1/L^2) and by the number of clusters
// used for the track fitting (1/sqrt(Ncl)). 
// Combining/refitting  the two segments together significantly better momentum resolution can be obtained.
//    sigma_{1/pt} ~ 8x10^-3  - defaul tracking (e.g only upper track)
//    sigma_{1/pt} ~ 8x10^-4  - combined tracking
//===========================================================================================
// 
// Interface/Implementation:
// This class provide the functionality to combined two ESD tracks and perform the trackfit using all related track points.
// Input: AliESDtrack's (or AliTPCseed's, expert use)
// Output: the trackparams from the combined fit are stored as data member of this class, and they can be accessed via getter or copy.
//
//===========================================================================================      
//
// Algorithm:
// 1. Read in two ESDtracks (see function CombineESDtracks())
//      TPCseeds are got from ESDfriends
//      Cuts are applied on the following quantities: #TPC cluster, relation between global-y coordinate of the two tracks, lever arm.
// 2. AliTPCCosmicUtils::CombinedFit is called and refit is performed using the TPCseeds from both tracks
// 3. The fitting quality is analyzed and stored and can be accessed via function GetStatus(). The detailed judgement of the quality is recorded as "enum CombineStatus"
//
//
//===========================================================================================
// Algorithm numerical debugging:
// The class can be used in the different debug/verbose level (see fDebugLevel) 
// Several intermediate variables can be stored in the trees, printout, or draw.
// Given functionality (dumping of variables to the tree) was also used for the tuning of the pair
// selection criterias, and for validation of the fit functionality. 
//
// ----- Debug output:
// for (debuglevel & 1)==1 && good fit, the following info saved:
/*
  (*fStreamer)<<"TrackProp"<<
  "Tup.="<<fTrackparUp<<          //AliExternalTrackParam at uppermost cluster obtained by upward propagation
  "Tlow.="<<fTrackparLow<<        //AliExternalTrackParam at lowermost cluster obtained by downward propagation
  "icup.="<<&fInnerClusterUp<<    //TVector3 position of the innermost cluster of the upper track
  "iclow.="<<&fInnerClusterLow<<
  "leverarm="<<fLeverArm<<
  "ncl="<<fFitNcls<<              //number of clusters used in successful propagation  
  "nmiss="<<fMissNcls<<           //number of clusters failed in propagation, should always be 0 in this case.
  "chi2="<<fPreChi2<<             //chi2/nfit  
  "momup="<<  momup <<            //momentum at uppermost cluster with upward propagation
  "momlow="<< momlow <<           //momentum at lowermost cluster with downward propagation
  "ptup="<<   ptup <<
  "ptlow="<<  ptlow <<
      "\n";
*/
// for (debuglevel & 2)==1, debug info in AliTPCCosmicUtils::FitKernel saved
//
//===========================================================================================
// Usage
/*
  fCombinedTrackfit = new AliTPCCosmicTrackfit(debuglevel, "anystring");

  //order not important; will be internally ordered (potinters modified due to &) such that track0 is the upper one
  //kfit = kTRUE: good fit, kFALSE: bad fit
  const Bool_t kfit = fCombinedTrackfit->CombineESDtracks(esdtrack0, esdtrack1);

  //status = 0 for good fit (i.e. kfit=kTRUE), non-0 for bad fit (i.e. kfit=kFALSE), see "enum CombineStatus" definition in header file
  const Int_t status = fCombinedTrackfit->GetStatus(); 
*/
//===========================================================================================
// Efficiency
//
// for 2011 Feb. cosmic data nch=2 events, the kfit and status look like:
/*
kfit,status (  0,   1):   16337 /  443901 =   3.680%          //kFailGetTPCseeds
kfit,status (  0,   2):    3692 /  443901 =   0.832%          //not both tracks have ncl > AliTPCCosmicUtils::fgkNclsMin
kfit,status (  0,   3):    6527 /  443901 =   1.470%          //clusters in two tracks should be clearly separated in y, i.e. lowest cluster of upper track higher than highest cluster of lower track; otherwise fail 
kfit,status (  0,   4):    7033 /  443901 =   1.584%          //fLeverArm<fgkCutLeverArm
kfit,status (  0,   6):    4398 /  443901 =   0.991%          //fail in propagation of at least one cluster
kfit,status (  0,   7):     508 /  443901 =   0.114%          //chi2/nfit > fgkMaxChi2
kfit,status (  0,   8):      52 /  443901 =   0.012%          //fail in geting impact parameters
kfit,status (  1,   0):  405354 /  443901 =  91.316%          //i.e. 91% of nch=2 events are successfully fitted.
*/
//
//===========================================================================================
// Resolution
//
// for muon momentum small than 20 GeV, energy loss in muon filter is visable when compaing fTrackparUp and fTrackparLow; energy loss estimated as 5 MeV/cm.
// particle traversing muon filter can be rejected by requiring "fInnerClusterUp.fZ > -40 && fInnerClusterLow.fZ > -40"
// momentum resolution is estimated by comparing the trackfit result by upward propagation through odd pad rows and that by downward propagation through even pad rows. Number of clusters used in this case is only half of that in normal usage.
// RMS of log10 p = 0.01 at 10 GeV/c, 0.1 at 100 GeV/c, 0.5 at 1 TeV/c.
// muon filter deteriorates momentum resolution by about +0.01 (absolute value).
//

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
#include "AliTPCCosmicTrackfit.h"

AliTPCCosmicTrackfit::AliTPCCosmicTrackfit(const Int_t dlev, const TString tag):
  fStreamer(0x0), fDebugLevel(dlev)
  , fSeedUp(0x0), fSeedLow(0x0), fTrackparUp(0x0), fTrackparLow(0x0), fIniTrackparUp(0x0), fIniTrackparLow(0x0)
  , fStatus(-999)
  , fLeverArm(-999)
  , fFitNcls(-999), fMissNcls(-999), fPreChi2(-999), fFitLeverArm(-999), fImpactD(-999), fImpactZ(-999)
{
  //
  //Constructor
  //
  fInnerClusterUp.SetXYZ(-999,-999,-999);
  fInnerClusterLow.SetXYZ(-999,-999,-999);

  if(fDebugLevel>0)
    fStreamer = new TTreeSRedirector(Form("CombinedTrackfit_%s.root", tag.Data()));

  SetRow(0,1);
  SetX(-1e10, 1e10);
}

AliTPCCosmicTrackfit::~AliTPCCosmicTrackfit()
{
  //
  //Destructor
  //

  delete fStreamer;
  
  delete fTrackparUp;
  delete fTrackparLow;
  delete fIniTrackparUp;
  delete fIniTrackparLow;
}

Bool_t AliTPCCosmicTrackfit::CombineESDtracks(AliESDtrack * &trk0, AliESDtrack *&trk1)
{
  //
  //Get TPCseeds from the 2 ESDtracks, swap TPCseeds and ESDTracks (if necessary) according to y (0:upper 1:lower), perform trackfit using TPCseeds
  //if fStatus==0, i.e. combine is successful, swap of the ESDtracks is kept since pointer *& is used
  //

  IniCombineESDtracks();

  if(!GetTPCseeds(trk0, trk1)){
    return kFALSE; 
  }

  CombineTPCseeds();

  if(fStatus == 0){
    if(fKswap){
      AliESDtrack * tmptrk = trk0;
      trk0 = trk1;
      trk1 = tmptrk;
    }
    return kTRUE;
  }
  else
    return kFALSE;
}

Bool_t AliTPCCosmicTrackfit::CombineTPCseedsFast(AliTPCseed * tpcseeds[], const AliExternalTrackParam * trkpars[])
{
  //
  //do combined track fit for given TPC seeds and initial trackpar, [0]: upper, [1]: lower
  //

  IniCombineESDtracks();

  fSeedUp = tpcseeds[0];
  fSeedLow = tpcseeds[1];

  fTrackparUp =  new AliExternalTrackParam(*(trkpars[0]));
  fTrackparLow = new AliExternalTrackParam(*(trkpars[1]));

  AliExternalTrackParam * trackPars[]={fTrackparUp, fTrackparLow};
  const AliTPCseed *seeds[]={fSeedUp, fSeedLow};
  TTreeSRedirector * debugstreamer = 0x0;
  if(fDebugLevel&2){
    debugstreamer = fStreamer;
  }

  AliTPCCosmicUtils::CombinedFit(trackPars, seeds, fRowStartShift, fRowStep, fXMin, fXMax, fFitNcls, fMissNcls, fPreChi2, fFitLeverArm, fImpactD, fImpactZ, debugstreamer);

  Update();

  if(fStatus == 0){
    return kTRUE;
  }
  else
    return kFALSE;
}

Bool_t AliTPCCosmicTrackfit::CombineTPCseeds(AliTPCseed * &seed0, AliTPCseed *&seed1)
{
  //
  //same as AliTPCCosmicTrackfit::CombineESDtracks, except that the seeds are passed in from outside, which can be still unordered
  //if fStatus==0, i.e. combine is successful, swap of the TPCseeds is kept since pointer *& is used
  //
  IniCombineESDtracks();

  fSeedUp  = seed0;
  fSeedLow = seed1;
  
  CombineTPCseeds();

  if(fStatus==0){
    if(fKswap){
      AliTPCseed * tmpseed = seed0;
      seed0 = seed1;
      seed1 = tmpseed;
    }
    return kTRUE;
  }
  else 
    return kFALSE;
}

void AliTPCCosmicTrackfit::Print() const
{
  //
  //print out variable values
  //
  printf("Status %2d NclsU %3d NclsD %3d ZinnerU %7.2f ZinnerD %7.2f LeverArm %7.2f\n", fStatus, fSeedUp->GetNumberOfClusters(), fSeedLow->GetNumberOfClusters(), fInnerClusterUp.Z(), fInnerClusterLow.Z(), fLeverArm);
}

/*
Double_t AliTPCCosmicTrackfit::ImpactParameter2D() const
{
  //
  //calculate the 2D-impactparameter from (0,0)
  //

  const TVector3 p0(0,0,0);
  const TVector3 l1(fInnerClusterUp.X(),  fInnerClusterUp.Y(),  0);
  const TVector3 l2(fInnerClusterLow.X(), fInnerClusterLow.Y(), 0);

  return AliTPCCosmicUtils::Point2LineDist(p0, l1, l2);
}

Double_t AliTPCCosmicTrackfit::ImpactParameter3D() const
{
  //
  //calculate the 3D-impactparameter from (0,0,0)
  //

  const TVector3 p0(0,0,0);
  const TVector3 l1(fInnerClusterUp);
  const TVector3 l2(fInnerClusterLow);

  return AliTPCCosmicUtils::Point2LineDist(p0, l1, l2);
}
*/

Double_t AliTPCCosmicTrackfit::MinPhi() const
{
  //
  //the smaller phi of the two tracks w.r.t. horizon
  //
  Double_t fsp[] = {TMath::Abs(TMath::Sin(fTrackparUp->Phi())), TMath::Abs(TMath::Sin(fTrackparLow->Phi()))};;
  return TMath::ASin(TMath::Min(fsp[0], fsp[1])) * TMath::RadToDeg();
}
//===================================================================================================
//===================================================================================================

void AliTPCCosmicTrackfit::IniCombineESDtracks()
{
  //
  //initialization, for reuse of the same AliTPCCosmicTrackfit instance
  //

  fSeedUp = 0x0;
  fSeedLow = 0x0;
  delete fTrackparUp;
  delete fTrackparLow;
  fTrackparUp = 0x0;
  fTrackparLow = 0x0;

  delete fIniTrackparUp;
  delete fIniTrackparLow;
  fIniTrackparUp = 0x0;
  fIniTrackparLow = 0x0;

  fStatus = 0;
  fKswap = kFALSE;
}

void AliTPCCosmicTrackfit::CombineTPCseeds()
{
  //
  //do combined trackfit using TPCseeds
  //

  if(
     !CheckNcls()
     || !AnaSeeds() // fSeedUp/Low swapped here! check by runTest.C 
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

  fIniTrackparUp = new AliExternalTrackParam(*fTrackparUp);
  fIniTrackparLow = new AliExternalTrackParam(*fTrackparLow);

  AliExternalTrackParam * trackPars[]={fTrackparUp, fTrackparLow};
  const AliTPCseed *seeds[]={fSeedUp, fSeedLow};
  TTreeSRedirector * debugstreamer = 0x0;
  if(fDebugLevel&2){
    debugstreamer = fStreamer;
  }

  AliTPCCosmicUtils::CombinedFit(trackPars, seeds, fRowStartShift, fRowStep, fXMin, fXMax, fFitNcls, fMissNcls, fPreChi2, fFitLeverArm, fImpactD, fImpactZ, debugstreamer);

  Update();

  return;
}

void AliTPCCosmicTrackfit::Update()
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

  if(fImpactD<0){
    fStatus = kFailImpact;
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

Bool_t AliTPCCosmicTrackfit::CheckLeverArm()
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

Bool_t AliTPCCosmicTrackfit::AnaSeeds()
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
    fKswap = kFALSE;

    kpass = kTRUE;

    combinedmax = singlemax[0];
    combinedmin = singlemin[1];
  }
  else if(singlemin[1].Y()> singlemax[0].Y()){
    fInnerClusterUp  = singlemin[1];
    fInnerClusterLow = singlemax[0];
  
    //have to be swapped
    fKswap = kTRUE;
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

Bool_t AliTPCCosmicTrackfit::CheckNcls()
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

Bool_t AliTPCCosmicTrackfit::GetTPCseeds(const AliESDtrack *trk0,  const AliESDtrack *trk1)
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

