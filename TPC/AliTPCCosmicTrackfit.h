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
//  Xianguo Lu 
//  lu@physi.uni-heidelberg.de
//  Xianguo.Lu@cern.ch

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

class AliTPCCosmicTrackfit
{
 public:
  AliTPCCosmicTrackfit(const Int_t dlev=0, const TString tag="test");
  
  ~AliTPCCosmicTrackfit();

  void SetRow(const Int_t shift, const Int_t step){ fRowStartShift = shift; fRowStep = step; }
  void SetX(const Double_t xmin, const Double_t xmax){ fXMin = xmin; fXMax = xmax; }

  Bool_t CombineESDtracks(AliESDtrack * &trk0, AliESDtrack *&trk1);
  Bool_t CombineTPCseedsFast(AliTPCseed * tpcseeds[], const AliExternalTrackParam * trkpars[]);
  Bool_t CombineTPCseeds(AliTPCseed * &seed0, AliTPCseed *&seed1);
  void Print() const ;
  
  //--------- getters ------------
  Bool_t IsSwap() const {return fKswap;}
  Int_t GetStatus()const{return fStatus;}
  Int_t GetFitNcls()const{return fFitNcls;}
  Int_t GetMissNcls()const{return fMissNcls;}
  Double_t GetChi2PerCluster()const{return fPreChi2;}
  Double_t GetFitLeverArm() const {return fFitLeverArm;}
  Double_t GetImpactD() const {return fImpactD;}
  Double_t GetImpactZ() const {return fImpactZ;}
  Double_t GetLeverArm()const {return fLeverArm;}
  TVector3 GetInnerClusterUp()const {return fInnerClusterUp;}
  TVector3 GetInnerClusterLow()const {return fInnerClusterLow;}
  /*
  Double_t ImpactParameter2D() const;
  Double_t ImpactParameter3D() const;
  */
  Double_t MinPhi()const;

  const AliExternalTrackParam * GetTrackParamUp() const {return fTrackparUp;}
  const AliExternalTrackParam * GetTrackParamLow() const {return fTrackparLow;}
  const AliExternalTrackParam * GetIniTrackParamUp() const {return fIniTrackparUp;}
  const AliExternalTrackParam * GetIniTrackParamLow() const {return fIniTrackparLow;}

  AliExternalTrackParam * CopyTrackParamUp() const {return new AliExternalTrackParam(*fTrackparUp);}
  AliExternalTrackParam * CopyTrackParamLow() const {return new AliExternalTrackParam(*fTrackparLow);}
  AliExternalTrackParam * CopyIniTrackParamUp() const {return new AliExternalTrackParam(*fIniTrackparUp);}
  AliExternalTrackParam * CopyIniTrackParamLow() const {return new AliExternalTrackParam(*fIniTrackparLow);}

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
    kFailChi2       =7,
    kFailImpact     =8
  };
  
  AliTPCCosmicTrackfit(const AliTPCCosmicTrackfit & p);
  AliTPCCosmicTrackfit & operator=(const AliTPCCosmicTrackfit & p);

  void IniCombineESDtracks();
  Bool_t GetTPCseeds(const AliESDtrack *trk0,  const AliESDtrack *trk1);
  Bool_t CheckNcls();
  Bool_t CheckLeverArm();
  Bool_t AnaSeeds();

  void CombineTPCseeds();
  void Update();

  TTreeSRedirector *fStreamer;     //!debug streamer
  Int_t fDebugLevel;                    //debug level

  AliTPCseed * fSeedUp;                         //TPC seed of upper track
  AliTPCseed * fSeedLow;                        //TPC seed of lower track
  AliExternalTrackParam * fTrackparUp;          //track param of upper track
  AliExternalTrackParam * fTrackparLow;         //track param of lower track
  AliExternalTrackParam * fIniTrackparUp;          //track param of upper track by MakeSeed
  AliExternalTrackParam * fIniTrackparLow;         //track param of lower track by MakeSeed

  Int_t fStatus;                               //status for CombineESDtracks/CombineTPCseeds: 0-successful, otherwise fail
  Bool_t fKswap;                               //true if should be / already be swapped

  TVector3 fInnerClusterUp;                    //xyz of the inner most TPC trackpoint of the Upper track
  TVector3 fInnerClusterLow;                   //xyz of the inner most TPC trackpoint of the Lower track
  Double_t fLeverArm;                          //transverse difference between upper most and lower most clusters

  Int_t fFitNcls;                              //number of TPC clusters successful in propagation (mean of the two propagation: upwards and downwards)
  Int_t fMissNcls;                             //number of TPC clusters fail in propagation (sum of the two propagation)
  Double_t fPreChi2;                           //Predicted chi2/nfit over the two propagation
  Double_t fFitLeverArm;                       //Lever arm calculated from clusters in fitkernel
  Double_t fImpactD;                           //2d impact parameter
  Double_t fImpactZ;                           //z of impact parameter

  Int_t fRowStartShift;                        //row start shift
  Int_t fRowStep;                              //row step
  Double_t fXMin;                              //cluster X min
  Double_t fXMax;                              //cluster X max
  static const Double_t fgkCutLeverArm = 350;  //minimum lever arm 350 ~ 250 * sqrt(2)
  static const Double_t fgkMaxChi2 = 10;       //max. chi2/ncls
};

#endif
