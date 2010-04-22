//$Id  $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Svein Lindal (slindal@fys.uio.no)                 *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTGlobalTrackMatcher.cxx
    @author Svein Lindal
    @date   
    @brief  The HLT class Matching Calorimeter clusters to TPC tracks
*/

#include "AliHLTGlobalTrackMatcher.h"


#if __GNUC__>= 3
using namespace std;
#endif

ClassImp(AliHLTGlobalTrackMatcher)


AliHLTGlobalTrackMatcher::AliHLTGlobalTrackMatcher() :
  fMaxZ(0),
  fMaxX(0),
  fMatchDistance(0),
  fRadius(0),
  fYSign(kFALSE)
{
  //Default constructor
}

AliHLTGlobalTrackMatcher::~AliHLTGlobalTrackMatcher()
{
  //Destructor
}

Int_t AliHLTGlobalTrackMatcher::AddTrackToCluster(Int_t tId, TArrayI* matchedTracksArray, Bool_t bestMatch, Int_t nMatches){
  //See header file for documentation
    
  matchedTracksArray->Set(matchedTracksArray->GetSize() + 1);
  if ( bestMatch )  {
    matchedTracksArray->AddAt(matchedTracksArray->At(0), matchedTracksArray->GetSize() - 1);
    matchedTracksArray->AddAt(tId, 0);
  } else {
    matchedTracksArray->AddAt(tId, matchedTracksArray->GetSize() - 1);
  }

  return nMatches;
}

Int_t AliHLTGlobalTrackMatcher::AddTrackToCluster(Int_t tId, Int_t* matchArray, bool bestMatch, Int_t nMatches ){
  //See header file for documentation
  
  //BALLE TODO: remove hardcoded 9
  if (nMatches > 9) {                                                   //BALLE this one tooo
    HLTError("The number of matching tracks (%d) exceeds the array size of %d", nMatches, 9);
    return 0;
  }

  if(bestMatch) {
    matchArray[nMatches] = matchArray[0];
    matchArray[0] = tId;
  } else  {
    matchArray[nMatches] = tId;
  }
  
  return nMatches;

};


Bool_t AliHLTGlobalTrackMatcher::IsTrackCloseToDetector(AliExternalTrackParam * track, Double_t bz, Double_t fMaxX, Bool_t ySign, Double_t fMaxZ, Double_t dRadius) {
  //See header file for documentation

  //Get track instersection with cylinder defined by detector radius
  Double_t trackPosition[3] = {0, 0, 0};
  if (! (track->GetXYZAt(dRadius, bz, trackPosition)) ) {
    return kFALSE;
  }

  //Positive y for EMCAL, negative for PHOS
  if(ySign) {
    if (trackPosition[1] < 0 ) {
      return kFALSE;
    }
    
  } else {
    if (trackPosition[1] > 0 ) {
      return kFALSE;
    }
  }
  
  
  if ( (TMath::Abs(trackPosition[2]) > fMaxZ) ) 
    return kFALSE;
  
  if (TMath::Abs(trackPosition[0]) > fMaxX )
    return kFALSE;


  return kTRUE;  
}
