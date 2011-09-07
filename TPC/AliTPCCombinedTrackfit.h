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
//
// ----- Usage:
/*
  fCombinedTrackfit = new AliTPCCombinedTrackfit(debuglevel, "anystring");

  //order not important; will be internally ordered (potinters modified due to &) such that track0 is the upper one
  //kfit = kTRUE: good fit, kFALSE: bad fit
  const Bool_t kfit = fCombinedTrackfit->CombineESDtracks(esdtrack0, esdtrack1);

  //status = 0 for good fit (i.e. kfit=kTRUE), non-0 for bad fit (i.e. kfit=kFALSE), see "enum CombineStatus" definition in header file
  const Int_t status = fCombinedTrackfit->GetStatus(); 

  //in Analysis Task write when terminate 
  fCombinedTrackfit->GetStreamer()->GetFile()->Write();
*/
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
// ----- Efficiency:
// for 2011 Feb. cosmic data nch=2 events, the kfit and status look like:
/*
kfit,status (  0,   1):   68939 / 2611959 =   2.639%          //kFailGetTPCseeds
kfit,status (  0,   2):   14886 / 2611959 =   0.570%          //not both tracks have ncl > AliTPCCosmicUtils::fgkNclsMin
kfit,status (  0,   3):   53185 / 2611959 =   2.036%          //clusters in two tracks should be clearly separated in y, i.e. lowest cluster of upper track higher than highest cluster of lower track; otherwise fail
kfit,status (  0,   4):   39841 / 2611959 =   1.525%          //fLeverArm<fgkCutLeverArm
kfit,status (  0,   6):   12933 / 2611959 =   0.495%          //fail in propagation of at least one cluster
kfit,status (  0,   7):   19994 / 2611959 =   0.765%          //chi2/nfit > fgkMaxChi2
kfit,status (  1,   0): 2402181 / 2611959 =  91.969%          //i.e. 92% of nch=2 events are successfully fitted.
*/
//
// ----- Resolution:
// for muon momentum small than 20 GeV, energy loss in muon filter is visable when compaing fTrackparUp and fTrackparLow; energy loss estimated as 5 MeV/cm.
// particle traversing muon filter can be rejected by requiring "fInnerClusterUp.fZ > -40 && fInnerClusterLow.fZ > -40"
// momentum resolution is estimated by comparing the trackfit result by upward propagation through odd pad rows and that by downward propagation through even pad rows. Number of clusters used in this case is only half of that in normal usage.
// RMS of log10 p = 0.01 at 10 GeV/c, 0.1 at 100 GeV/c, 0.5 at 1 TeV/c.
// muon filter deteriorates momentum resolution by about +0.01 (absolute value).
//
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//
/*
//in [cm]
const Double_t _TPCZMIN = -250;
const Double_t _TPCZMAX =  250;
const Double_t _TPCINR  =  84.8;
const Double_t _TPCOUR  =  246.6;
*/

#ifndef ALITPCCOMBINEDTRACKFIT_H
#define ALITPCCOMBINEDTRACKFIT_H

class TTreeSRedirector;

class AliTPCCombinedTrackfit
{
 public:
  AliTPCCombinedTrackfit(const Int_t dlev=0, const TString tag="test");
  
  ~AliTPCCombinedTrackfit();

  Bool_t CombineESDtracks(AliESDtrack * &trk0, AliESDtrack *&trk1);
  Bool_t CombineTPCseeds(AliTPCseed * &seed0, AliTPCseed *&seed1);
  void Print() const ;
  //--------- getters ------------

  Int_t GetStatus()const{return fStatus;}
  Int_t GetFitNcls()const{return fFitNcls;}
  Int_t GetMissNcls()const{return fMissNcls;}
  Double_t GetChi2PerCluster()const{return fPreChi2;}
  Double_t GetLeverArm()const {return fLeverArm;}
  TVector3 GetInnerClusterUp()const {return fInnerClusterUp;}
  TVector3 GetInnerClusterLow()const {return fInnerClusterLow;}
  Double_t ImpactParameter() const;
  Double_t MinPhi()const;

  AliExternalTrackParam * GetTrackParamUp() const {return fTrackparUp;}
  AliExternalTrackParam * GetTrackParamLow() const {return fTrackparLow;}
  AliTPCseed * GetTPCseedUp()   const {return fSeedUp;}
  AliTPCseed * GetTPCseedLow() const {return fSeedLow;}

  TTreeSRedirector * GetStreamer() const {return fStreamer;}

 private:
  enum CombineStatus{
    kFailGetTPCseeds=1,
    kFailNclsMin    =2,
    kFailSwapSeeds  =3,
    kFailLeverArm   =4,
    kFailMakeSeed   =5,
    kFailPropagation=6,
    kFailChi2       =7
  };
  
  AliTPCCombinedTrackfit(const AliTPCCombinedTrackfit & p);
  AliTPCCombinedTrackfit & operator=(const AliTPCCombinedTrackfit & p);

  void IniCombineESDtracks();
  Bool_t GetTPCseeds(const AliESDtrack *trk0,  const AliESDtrack *trk1);
  Bool_t CheckNcls();
  Bool_t CheckLeverArm();
  Bool_t AnaSeeds(Bool_t &kswap);

  void CombineTPCseeds(Bool_t &kswap);
  void Update();

  TTreeSRedirector *fStreamer;     //!debug streamer
  Int_t fDebugLevel;                    //debug level

  AliTPCseed * fSeedUp;                         //TPC seed of upper track
  AliTPCseed * fSeedLow;                        //TPC seed of lower track
  AliExternalTrackParam * fTrackparUp;          //track param of upper track
  AliExternalTrackParam * fTrackparLow;         //track param of lower track

  Int_t fStatus;                               //status for CombineESDtracks/CombineTPCseeds: 0-successful, otherwise fail
  
  TVector3 fInnerClusterUp;                    //xyz of the inner most TPC trackpoint of the Upper track
  TVector3 fInnerClusterLow;                   //xyz of the inner most TPC trackpoint of the Lower track
  Double_t fLeverArm;                          //transverse difference between upper most and lower most clusters

  Int_t fFitNcls;                              //number of TPC clusters successful in propagation (mean of the two propagation: upwards and downwards)
  Int_t fMissNcls;                             //number of TPC clusters fail in propagation (sum of the two propagation)
  Double_t fPreChi2;                           //Predicted chi2/nfit over the two propagation

  static const Double_t fgkCutLeverArm = 350;  //minimum lever arm 350 ~ 250 * sqrt(2)
  static const Double_t fgkMaxChi2 = 10;       //max. chi2/ncls
};

#endif
