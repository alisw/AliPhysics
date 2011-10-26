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
//comment
//comment
//comment
//comment
//comment
//
//
//  Xianguo Lu 
//  lu@physi.uni-heidelberg.de
//  Xianguo.Lu@cern.ch
//


#ifndef ALIESDCOSMICTRACK_H
#define ALIESDCOSMICTRACK_H

//#include "AliExternalTrackParam.h"

class AliESDCosmicTrack: public AliExternalTrackParam
{
 public:
  AliESDCosmicTrack();
  AliESDCosmicTrack(const Int_t idUp, const Int_t idLow, const AliExternalTrackParam * trkparUp, const AliExternalTrackParam * trkparLow, const AliExternalTrackParam * parx0Up, const AliExternalTrackParam * parx0Low, const Int_t ncls, const Double_t la, const Double_t chi2, const Double_t vd, const Double_t vz, const Bool_t isreuse, const Double_t findable);
  AliESDCosmicTrack(const AliESDCosmicTrack & costrk);
  AliESDCosmicTrack operator=(const AliESDCosmicTrack & costrk);

  virtual ~AliESDCosmicTrack(){}

  Int_t GetESDUpperTrackIndex() const {return fESDtrackIndex[0];}
  Int_t GetESDLowerTrackIndex() const {return fESDtrackIndex[1];}
  const AliExternalTrackParam * GetLowerPartner() const {return &fLowerTrackParam;}
  const AliExternalTrackParam * GetESDUpperTrackParamAt0() const {return &fX0UpperTrackParam;}
  const AliExternalTrackParam * GetESDLowerTrackParamAt0() const {return &fX0LowerTrackParam;}

  Int_t GetNCluster() const {return fNCluster;}
  Double_t GetLeverArm() const {return fLeverArm;}
  Double_t GetChi2PerCluster() const {return fChi2PerCluster;}
  Double_t GetImpactD() const {return fImpactD;}
  Double_t GetImpactZ() const {return fImpactZ;}
  Bool_t IsReuse() const{return fIsReuse;}
  Double_t GetMinFindableRatio()const{return fFindableRatio;}

 private:
  Int_t fESDtrackIndex[2];                                    //[0]= ESD track for this object; [1]= ESD track for LowerTrackParma

  AliExternalTrackParam fLowerTrackParam;                     //trackparameter estimated at lower-outer part of TPC
  AliExternalTrackParam fX0UpperTrackParam;                   //ESD upper trackparam at X=0
  AliExternalTrackParam fX0LowerTrackParam;                   //ESD lower trackparam at X=0


  Int_t fNCluster;                                            //number of cls used in fit
  Double_t fLeverArm;                                         //lever arm
  Double_t fChi2PerCluster;                                   //chi2/ncls
  Double_t fImpactD;                                          //2d impact parameter
  Double_t fImpactZ;                                          //z of impact parameter
  Bool_t fIsReuse;                                            //true if one of the track from the pair already used in previous pair
  Double_t fFindableRatio;                                    //min of TPC ncls/nfindablecls of the two tracks

  ClassDef(AliESDCosmicTrack, 1); //0: only data members derived from AliExternalTrackParam are saved in tree; 1: all are saved!!
};

#endif


