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
AliMUONGlobalTrigger::AliMUONGlobalTrigger(const AliMUONGlobalTrigger& MUONGlobalTrig):TObject(MUONGlobalTrig)
{
// copy constructor
  fSinglePlusLpt  = MUONGlobalTrig.fSinglePlusLpt;
  fSinglePlusHpt  = MUONGlobalTrig.fSinglePlusHpt;
  fSinglePlusApt  = MUONGlobalTrig.fSinglePlusApt;
  
  fSingleMinusLpt = MUONGlobalTrig.fSingleMinusLpt;
  fSingleMinusHpt = MUONGlobalTrig.fSingleMinusHpt;
  fSingleMinusApt = MUONGlobalTrig.fSingleMinusApt;
  
  fSingleUndefLpt = MUONGlobalTrig.fSingleUndefLpt;
  fSingleUndefHpt = MUONGlobalTrig.fSingleUndefHpt;
  fSingleUndefApt = MUONGlobalTrig.fSingleUndefApt;
  
  fPairUnlikeLpt  = MUONGlobalTrig.fPairUnlikeLpt;
  fPairUnlikeHpt  = MUONGlobalTrig.fPairUnlikeHpt;
  fPairUnlikeApt  = MUONGlobalTrig.fPairUnlikeApt;
  
  fPairLikeLpt    = MUONGlobalTrig.fPairLikeLpt;
  fPairLikeHpt    = MUONGlobalTrig.fPairLikeHpt;
  fPairLikeApt    = MUONGlobalTrig.fPairLikeApt;
}

//----------------------------------------------------------------------
AliMUONGlobalTrigger& AliMUONGlobalTrigger::operator=(const AliMUONGlobalTrigger& MUONGlobalTrig)
{
// equal operator (useful for non-pointer member in TClonesArray)
  if (this == &MUONGlobalTrig)
    return *this;

  fSinglePlusLpt  = MUONGlobalTrig.fSinglePlusLpt;
  fSinglePlusHpt  = MUONGlobalTrig.fSinglePlusHpt;
  fSinglePlusApt  = MUONGlobalTrig.fSinglePlusApt;
  
  fSingleMinusLpt = MUONGlobalTrig.fSingleMinusLpt;
  fSingleMinusHpt = MUONGlobalTrig.fSingleMinusHpt;
  fSingleMinusApt = MUONGlobalTrig.fSingleMinusApt;
  
  fSingleUndefLpt = MUONGlobalTrig.fSingleUndefLpt;
  fSingleUndefHpt = MUONGlobalTrig.fSingleUndefHpt;
  fSingleUndefApt = MUONGlobalTrig.fSingleUndefApt;
  
  fPairUnlikeLpt  = MUONGlobalTrig.fPairUnlikeLpt;
  fPairUnlikeHpt  = MUONGlobalTrig.fPairUnlikeHpt;
  fPairUnlikeApt  = MUONGlobalTrig.fPairUnlikeApt;
  
  fPairLikeLpt    = MUONGlobalTrig.fPairLikeLpt;
  fPairLikeHpt    = MUONGlobalTrig.fPairLikeHpt;
  fPairLikeApt    = MUONGlobalTrig.fPairLikeApt;

  return *this;
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

//----------------------------------------------------------------------
//--- methods which return member data related info
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::SinglePlusLpt(){
// returns Number of Single Plus Low pt 
  return fSinglePlusLpt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::SinglePlusHpt(){  
// returns Number of Single Plus High pt 
  return fSinglePlusHpt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::SinglePlusApt(){  
// returns Number of Single Plus All pt 
  return fSinglePlusApt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::SingleMinusLpt(){ 
// returns Number of Single Minus Low pt
  return fSingleMinusLpt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::SingleMinusHpt(){
// returns Number of Single Minus High pt 
  return fSingleMinusHpt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::SingleMinusApt(){
// returns Number of Single Minus All pt
  return fSingleMinusApt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::SingleUndefLpt(){
// returns Number of Single Undefined Low pt
  return fSingleUndefLpt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::SingleUndefHpt(){ 
// returns Number of Single Undefined High pt 
  return fSingleUndefHpt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::SingleUndefApt(){
// returns Number of Single Undefined All pt
  return fSingleUndefApt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::PairUnlikeLpt(){  
// returns Number of Unlike sign pair Low pt
  return fPairUnlikeLpt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::PairUnlikeHpt(){
// returns Number of Unlike sign pair High pt
  return fPairUnlikeHpt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::PairUnlikeApt(){
// returns Number of Unlike sign pair All pt
  return fPairUnlikeApt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::PairLikeLpt(){
// returns Number of Like sign pair Low pt
  return fPairLikeLpt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::PairLikeHpt(){
// returns Number of Like sign pair High pt
  return fPairLikeHpt;
}
//----------------------------------------------------------------------
Int_t AliMUONGlobalTrigger::PairLikeApt(){
// returns Number of Like sign pair All pt
  return fPairLikeApt;
}







