//$Id$
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
  fPhosMaxZ(0),
  fPhosMaxX(0),
  fEmcalMaxZ(0),
  fEmcalMaxX(0),
  fMatchDistance(0),
  fMatchDistanceEMCal(0),
  fPhosRadius(460),
  fEmcalRadius(448),
  fStep(100.),
  fMass(0.139)
{
  //Default constructor
  DoInit();
}

//_____________________________________________________________________________
AliHLTGlobalTrackMatcher::~AliHLTGlobalTrackMatcher()
{
  //Destructor

}

void AliHLTGlobalTrackMatcher::DoInit( ) {
  //See header file for documentation
  //BALLE TODO: Change hardcoded values to something that is initialised through command line or something!!!


  fMatchDistance = 40*40;
  fMatchDistanceEMCal = 0.1; // EMCal EtaxPhi cut 

  fPhosMaxX = 355 + TMath::Sqrt(fMatchDistance) + 30;
  fPhosMaxZ = 64.+ TMath::Sqrt(fMatchDistance) + 30;

  fEmcalMaxZ = 350 + TMath::Sqrt(fMatchDistance) + 30;
  fEmcalMaxX = 3000;

  fStep=100.;// Step for EMCAL extrapolation
  fMass=0.139;// Mass for EMCAL extrapolation hipothesis

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

  //  HLTInfo("Adding track %d to cluster with %d previous matches", tId, nMatches);
  
  //BALLE TODO: remove hardcoded 9
  if (nMatches > 9) {                                                   //BALLE this on tooo
    HLTDebug("The number of matching tracks (%d) exceeds the array size of %d", nMatches, 9);
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
  
  //Positive y for EMCAL, negative for PHOS
  if(ySign) {
    //EMCAL
    if	(track->Pt()<1.)	return kFALSE;
    if	(track->Eta()>.8||track->Eta()<-.8)	return kFALSE;
    if	(track->Phi()>4.||track->Phi()<1.)	return kFALSE;   
  } else {
    //PHOS
      //Get track instersection with cylinder defined by detector radius
    Double_t trackPosition[3] = {0, 0, 0};
    if (! (track->GetXYZAt(dRadius, bz, trackPosition)) ) {
      return kFALSE;
    }
    if (trackPosition[1] > 0 ) 
      return kFALSE;
    if ( (TMath::Abs(trackPosition[2]) > fMaxZ) ) 
      return kFALSE;
    if (TMath::Abs(trackPosition[0]) > fMaxX )
      return kFALSE;
  }
  return kTRUE;  
}
