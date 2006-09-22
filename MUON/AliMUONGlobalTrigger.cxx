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
 

/* $Id$ */


#include "AliMUONGlobalTrigger.h"
#include <assert.h>
#include "AliLog.h"
#include "AliMUONLocalStruct.h"

/// \class AliMUONGlobalTrigger
/// Global Trigger algorithm data output.
/// Built from Local and Regional algorithms.                          \n 
/// Update for copy & assigment operator,
/// add SetGlobalPattern and GetGlobalPattern method for rawdata 
/// (Ch. Finck)
/// \author Ph. Crochet

/// \cond CLASSIMP
ClassImp(AliMUONGlobalTrigger)
/// \endcond

//----------------------------------------------------------------------
AliMUONGlobalTrigger::AliMUONGlobalTrigger()
  : TObject(),
    fSinglePlusLpt(0),
    fSinglePlusHpt(0),
    fSinglePlusApt(0),
  
    fSingleMinusLpt(0),
    fSingleMinusHpt(0),
    fSingleMinusApt(0),
  
    fSingleUndefLpt(0),
    fSingleUndefHpt(0),
    fSingleUndefApt(0),
    
    fPairUnlikeLpt(0),
    fPairUnlikeHpt(0),
    fPairUnlikeApt(0),
    
    fPairLikeLpt(0),
    fPairLikeHpt(0),
    fPairLikeApt(0)
{ 
  /// Default constructor 
}

//----------------------------------------------------------------------
AliMUONGlobalTrigger::AliMUONGlobalTrigger(const AliMUONGlobalTrigger& theMUONGlobalTrig)
  : TObject(theMUONGlobalTrig),
  
    fSinglePlusLpt(theMUONGlobalTrig.fSinglePlusLpt),
    fSinglePlusHpt(theMUONGlobalTrig.fSinglePlusHpt),
    fSinglePlusApt(theMUONGlobalTrig.fSinglePlusApt),
    
    fSingleMinusLpt(theMUONGlobalTrig.fSingleMinusLpt),
    fSingleMinusHpt(theMUONGlobalTrig.fSingleMinusHpt),
    fSingleMinusApt(theMUONGlobalTrig.fSingleMinusApt),
    
    fSingleUndefLpt(theMUONGlobalTrig.fSingleUndefLpt),
    fSingleUndefHpt(theMUONGlobalTrig.fSingleUndefHpt),
    fSingleUndefApt(theMUONGlobalTrig.fSingleUndefApt),
    
    fPairUnlikeLpt(theMUONGlobalTrig.fPairUnlikeLpt),
    fPairUnlikeHpt(theMUONGlobalTrig.fPairUnlikeHpt),
    fPairUnlikeApt(theMUONGlobalTrig.fPairUnlikeApt),
    
    fPairLikeLpt(theMUONGlobalTrig.fPairLikeLpt),
    fPairLikeHpt(theMUONGlobalTrig.fPairLikeHpt),
    fPairLikeApt(theMUONGlobalTrig.fPairLikeApt)
{
  /// Copy constructor
}

//----------------------------------------------------------------------
AliMUONGlobalTrigger::~AliMUONGlobalTrigger()
{
  /// Destructor
}

//----------------------------------------------------------------------
AliMUONGlobalTrigger& AliMUONGlobalTrigger::operator=(const AliMUONGlobalTrigger& theMUONGlobalTrig)
{
  /// Assignement operator;
  /// equal operator (useful for non-pointer member in TClonesArray)

  if (this == &theMUONGlobalTrig)
    return *this;
    
  // base class assignement
  TObject::operator=(theMUONGlobalTrig);

  fSinglePlusLpt  = theMUONGlobalTrig.fSinglePlusLpt;
  fSinglePlusHpt  = theMUONGlobalTrig.fSinglePlusHpt;
  fSinglePlusApt  = theMUONGlobalTrig.fSinglePlusApt;
  
  fSingleMinusLpt = theMUONGlobalTrig.fSingleMinusLpt;
  fSingleMinusHpt = theMUONGlobalTrig.fSingleMinusHpt;
  fSingleMinusApt = theMUONGlobalTrig.fSingleMinusApt;
  
  fSingleUndefLpt = theMUONGlobalTrig.fSingleUndefLpt;
  fSingleUndefHpt = theMUONGlobalTrig.fSingleUndefHpt;
  fSingleUndefApt = theMUONGlobalTrig.fSingleUndefApt;
  
  fPairUnlikeLpt  = theMUONGlobalTrig.fPairUnlikeLpt;
  fPairUnlikeHpt  = theMUONGlobalTrig.fPairUnlikeHpt;
  fPairUnlikeApt  = theMUONGlobalTrig.fPairUnlikeApt;
  
  fPairLikeLpt    = theMUONGlobalTrig.fPairLikeLpt;
  fPairLikeHpt    = theMUONGlobalTrig.fPairLikeHpt;
  fPairLikeApt    = theMUONGlobalTrig.fPairLikeApt;

  return *this;
}

//----------------------------------------------------------------------
AliMUONGlobalTrigger::AliMUONGlobalTrigger(Int_t *singlePlus, 
					   Int_t *singleMinus,
					   Int_t *singleUndef,
					   Int_t *pairUnlike, Int_t *pairLike)
  : TObject(),
  
    fSinglePlusLpt(singlePlus[0]),
    fSinglePlusHpt(singlePlus[1]),
    fSinglePlusApt(singlePlus[2]),

    fSingleMinusLpt(singleMinus[0]),
    fSingleMinusHpt(singleMinus[1]),
    fSingleMinusApt(singleMinus[2]),

    fSingleUndefLpt(singleUndef[0]),
    fSingleUndefHpt(singleUndef[1]),
    fSingleUndefApt(singleUndef[2]),

    fPairUnlikeLpt(pairUnlike[0]),
    fPairUnlikeHpt(pairUnlike[1]),
    fPairUnlikeApt(pairUnlike[2]),

    fPairLikeLpt(pairLike[0]),    
    fPairLikeHpt(pairLike[1]),  
    fPairLikeApt(pairLike[2])
  
{
  /// Set the Global Trigger object
}

//-----------------------------------------------------------
void AliMUONGlobalTrigger::SetGlobalPattern(Int_t gloTrigPat)
{
  /// Set class member from global pattern
  /// coming from rawdata
  /// [Hpt, Lpt] with [+, -, LS, US]

  fSinglePlusLpt = (gloTrigPat     ) & 0x1;
  fSinglePlusHpt = (gloTrigPat >> 1) & 0x1; 
  fSinglePlusApt = (gloTrigPat >> 2) & 0x1;

  fSingleMinusLpt = (gloTrigPat >> 3) & 0x1;
  fSingleMinusHpt = (gloTrigPat >> 4) & 0x1;
  fSingleMinusApt = (gloTrigPat >> 5) & 0x1; 

  fSingleUndefLpt = (gloTrigPat >> 6) & 0x1;
  fSingleUndefHpt = (gloTrigPat >> 7) & 0x1;
  fSingleUndefApt = (gloTrigPat >> 8) & 0x1;

  fPairUnlikeLpt = (gloTrigPat >> 9) & 0x1;
  fPairUnlikeHpt = (gloTrigPat >> 10) & 0x1;
  fPairUnlikeApt = (gloTrigPat >> 11) & 0x1;

  fPairLikeLpt   = (gloTrigPat >> 12) & 0x1;
  fPairLikeHpt   = (gloTrigPat >> 13) & 0x1;
  fPairLikeApt   = (gloTrigPat >> 14) & 0x1;

}
//-----------------------------------------------------------
void AliMUONGlobalTrigger::SetGlobalPattern(UShort_t globalResponse)
{
  /// Set class member from global response
  /// coming from trigger electronics
  /// should be unformized with rawdata (->oct 06)
  /// [Hpt, Lpt, Apt] with [+, -, US, LS]
  fSinglePlusLpt = ((globalResponse & 0xC0)  >>  6) == 2;
  fSinglePlusHpt = ((globalResponse & 0xC00) >> 10) == 2;
  fSinglePlusApt = ((globalResponse & 0xC)   >>  2) == 2;

  fSingleMinusLpt = ((globalResponse & 0xC0)  >>  6) == 1;
  fSingleMinusHpt = ((globalResponse & 0xC00) >> 10) == 1;
  fSingleMinusApt = ((globalResponse & 0xC)   >>  2) == 1;

  fSingleUndefLpt = ((globalResponse & 0xC0)  >>  6) == 3;
  fSingleUndefHpt = ((globalResponse & 0xC00) >> 10) == 3;
  fSingleUndefApt = ((globalResponse & 0xC)   >>  2) == 3;

  fPairUnlikeLpt = (globalResponse & 0x10)  >> 4;
  fPairUnlikeHpt = (globalResponse & 0x100) >> 8;
  fPairUnlikeApt = (globalResponse & 0x1);
  
  fPairLikeLpt = (globalResponse & 0x20)  >> 5;
  fPairLikeHpt = (globalResponse & 0x200) >> 9;
  fPairLikeApt = (globalResponse & 0x2)   >> 1;
  
}
//-----------------------------------------------------------
void AliMUONGlobalTrigger::SetFromGlobalResponse(UChar_t globalResponse)
{
  /// Set class members from global response
  /// coming from rawdata
  /// [US:2, LS:2, Single:2] with [Hpt, Lpt]
  /// remove Apt

  // don't have the information anymore of the sign
  fSinglePlusLpt = fSingleMinusLpt = globalResponse & 0x1;
  fSinglePlusHpt = fSingleMinusHpt = (globalResponse >> 1) & 0x1;

  fPairUnlikeLpt = (globalResponse >> 4)  & 0x1;
  fPairUnlikeHpt = (globalResponse >> 5)  & 0x1;
  
  fPairLikeLpt = (globalResponse >> 2)  & 0x1;
  fPairLikeHpt = (globalResponse >> 3)  & 0x1;
  
}
//-----------------------------------------------------------
Int_t AliMUONGlobalTrigger::GetGlobalPattern() const
{
  /// Global trigger pattern calculation
  /// from class member values

  Int_t gloTrigPat = 0;

  if (SinglePlusLpt())  gloTrigPat|= 0x1;
  if (SinglePlusHpt())  gloTrigPat|= 0x2;
  if (SinglePlusApt())  gloTrigPat|= 0x4;
 
  if (SingleMinusLpt()) gloTrigPat|= 0x8;
  if (SingleMinusHpt()) gloTrigPat|= 0x10;
  if (SingleMinusApt()) gloTrigPat|= 0x20;
 
  if (SingleUndefLpt()) gloTrigPat|= 0x40;
  if (SingleUndefHpt()) gloTrigPat|= 0x80;
  if (SingleUndefApt()) gloTrigPat|= 0x100;
 
  if (PairUnlikeLpt())  gloTrigPat|= 0x200;
  if (PairUnlikeHpt())  gloTrigPat|= 0x400;
  if (PairUnlikeApt())  gloTrigPat|= 0x800;

  if (PairLikeLpt())    gloTrigPat|= 0x1000;
  if (PairLikeHpt())    gloTrigPat|= 0x2000;
  if (PairLikeApt())    gloTrigPat|= 0x4000;

  return gloTrigPat;

}


//-----------------------------------------------------------
UChar_t AliMUONGlobalTrigger::GetGlobalResponse() const
{
  /// Global trigger response
  /// from class member values
  /// [US:2, LS:2, Single:2] with [Hpt, Lpt]
  /// remove Apt

  UChar_t response = 0;
  UChar_t respUS  = 0;
  UChar_t respLS  = 0;
  UChar_t respS  = 0;

  if (SinglePlusLpt() || SingleMinusLpt())  respS |= 0x1;
  if (SinglePlusHpt() || SingleMinusHpt())  respS |= 0x2;

  if (PairLikeLpt())    respLS |= 0x1;
  if (PairLikeHpt())    respLS |= 0x1;
  respLS <<= 2;

  if (PairUnlikeLpt())  respUS |= 0x1;
  if (PairUnlikeHpt())  respUS |= 0x2;
  respUS <<= 4;

  response = respUS | respLS | respS;

  return response;

}
//----------------------------------------------------------------------
void AliMUONGlobalTrigger::Print(Option_t* opt) const
{
  //
  // Printing Global Trigger information
  //
  TString sopt(opt);
  sopt.ToUpper();
  if ( sopt.Contains("FULL") ) { 

      printf("===================================================\n");
      printf(" Global Trigger output       Low pt  High pt   All\n");
      printf(" number of Single Plus      :\t");
      printf("%i\t%i\t%i\t",SinglePlusLpt(),SinglePlusHpt(),SinglePlusApt());
      printf("\n");
      
      printf(" number of Single Minus     :\t");
      printf("%i\t%i\t%i\t",SingleMinusLpt(),SingleMinusHpt(),SingleMinusApt());
      printf("\n");
      
      printf(" number of Single Undefined :\t"); 
      printf("%i\t%i\t%i\t",SingleUndefLpt(),SingleUndefHpt(),SingleUndefApt());
      printf("\n");
      
      printf(" number of UnlikeSign pair  :\t"); 
      printf("%i\t%i\t%i\t",PairUnlikeLpt(),PairUnlikeHpt(),PairUnlikeApt());
      printf("\n");
      
      printf(" number of LikeSign pair    :\t");  
      printf("%i\t%i\t%i\t",PairLikeLpt(),PairLikeHpt(),PairLikeApt());
      printf("\n");
      
      printf("===================================================\n");
  }  
}


