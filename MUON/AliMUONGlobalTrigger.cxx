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

/*

*/

#include "AliMUONGlobalTrigger.h"

ClassImp(AliMUONGlobalTrigger);
//----------------------------------------------------------------------
AliMUONGlobalTrigger::AliMUONGlobalTrigger()
{
// constructor 
  fSinglePlusLpt = 0;
  fSinglePlusHpt = 0;
  fSinglePlusApt = 0;
  
  fSingleMinusLpt = 0;
  fSingleMinusHpt = 0;
  fSingleMinusApt = 0;
  
  fSingleUndefLpt = 0;
  fSingleUndefHpt = 0;
  fSingleUndefApt = 0;
  
  fPairUnlikeLpt = 0;
  fPairUnlikeHpt = 0;
  fPairUnlikeApt = 0;
  
  fPairLikeLpt   = 0;
  fPairLikeHpt   = 0;
  fPairLikeApt   = 0;
}
//----------------------------------------------------------------------
AliMUONGlobalTrigger::AliMUONGlobalTrigger(Int_t *singlePlus, 
					   Int_t *singleMinus,
					   Int_t *singleUndef,
					   Int_t *pairUnlike, Int_t *pairLike)
{
// Set the Global Trigger object
  fSinglePlusLpt = singlePlus[0];
  fSinglePlusHpt = singlePlus[1];
  fSinglePlusApt = singlePlus[2];

  fSingleMinusLpt = singleMinus[0];
  fSingleMinusHpt = singleMinus[1];
  fSingleMinusApt = singleMinus[2];

  fSingleUndefLpt = singleUndef[0];
  fSingleUndefHpt = singleUndef[1];
  fSingleUndefApt = singleUndef[2];

  fPairUnlikeLpt = pairUnlike[0];
  fPairUnlikeHpt = pairUnlike[1];
  fPairUnlikeApt = pairUnlike[2];

  fPairLikeLpt   = pairLike[0];  
  fPairLikeHpt   = pairLike[1];  
  fPairLikeApt   = pairLike[2];  
}
