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
// Input: ESDevent
// Functionality: find ESDtrack pairs according to some criteria to form one cosmic ray; AliTPCCosmicTrackfit then performs the combined track fit for the pair
//
//  Xianguo Lu 
//  lu@physi.uni-heidelberg.de
//  Xianguo.Lu@cern.ch
//

#ifndef ALICOSMICTRACKER_H
#define ALICOSMICTRACKER_H

class TTreeSRedirector;
class AliTPCCosmicTrackfit;

class AliCosmicTracker
{
 public:

  typedef Bool_t (*CutFunc)(AliESDtrack *trk);

  AliCosmicTracker(const Int_t dlev=0, const TString tag="test"); 
  ~AliCosmicTracker();

  void SetESDEvent(AliESDEvent *esd);
  Int_t Process(const TString tag="", const Bool_t kprint=kFALSE);
  TClonesArray * GetTrackStack() const {return fTrackStack;}

  TTreeSRedirector * GetStreamer() const {return fStreamer;}

  void SetCutPull(     const Int_t ii, const Double_t cut){ fCutPull[ii] = cut;}
  void SetCutDelta(    const Int_t ii, const Double_t cut){ fCutDelta[ii] = cut;}

  void SetCutdPhi(const Double_t cut){fCutdPhi = cut;}
  void SetCutdTheta(const Double_t cut){fCutdTheta = cut;}

  void SetUserESDtrackCut(CutFunc func){fUserCut = func;}

 private:
  AliCosmicTracker(const AliCosmicTracker &p);
  AliCosmicTracker & operator=(const AliCosmicTracker &p);

  Bool_t ESDtrackCut(AliESDtrack * trk, Double_t &findabler);

  Bool_t IsPair(AliESDtrack* trk0, AliESDtrack*trk1);
  void WriteStreamer(Int_t ntrk);

  CutFunc fUserCut;                               //user ESDtrack function
  TTreeSRedirector *fStreamer;                   //debug streamer
  Int_t fDebugLevel;                             //debug level

  AliESDEvent *fESDEvent;                        //esd
  AliTPCCosmicTrackfit *fCombinedTrackfit;     //combinedtrackfit
  TClonesArray *fTrackStack;                     //storing cosmic ray

  AliExternalTrackParam fCombTrackUp;            //upper track param from combined track fit
  AliExternalTrackParam fCombTrackLow;          //lower track param from combined track fit
  AliExternalTrackParam fTrack0;                //upper track param estimated at x=0 from IsPair
  AliExternalTrackParam fTrack1;                //lower track param estimated at x=0 from IsPair

  Double_t fCombRefPt;                           //pt of the ray at lower-outer TPC by combined fit
  Double_t fChi2PerCluster;                      //chi2/ncls of the fit
  Double_t fImpactD;                                  //2d impact parameter
  Double_t fImpactZ;                                   //z of impact parameter
  Bool_t fIsReuse;                                     //if one the the tracks in the pair is already used in previous pair
  Double_t fFindableRatio;                      //min of TPC ncls/nfindablecls of the two tracks

  Double_t fdPhi;                                //phi0-phi1-pi of the EDS tracks
  Double_t fCutdPhi;                           //cut

  Double_t fdTheta;                              //theta0+theta1-pi of the ESD tracks
  Double_t fCutdTheta;                         //cut

  Double_t fPull[5];                             //pull of the two ESD tracks at x=0
  Double_t fCutPull[5];                            //pull

  Double_t fDelta[5];                             //delta of the two ESD tracks at x=0
  Double_t fCutDelta[5];                          //delta cut

  Int_t fErrFlagIsPair;                            //error status in IsPair

  static const Double_t fgkCutFindable = 0.5;            //cut on findable ratio of TPC cluster
};

#endif


