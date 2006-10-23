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


#include "AliMUONRegionalTrigger.h"
#include <assert.h>
#include "AliLog.h"

// -----------------------------
// Class AliMUONRegionalTrigger
// -----------------------------
// Regional Trigger algorithm data outputs
// Author:  Ch. Finck

/// \cond CLASSIMP
ClassImp(AliMUONRegionalTrigger)
/// \endcond

//----------------------------------------------------------------------
AliMUONRegionalTrigger::AliMUONRegionalTrigger()
  : TObject(), 
    fId(0),
    fLocalMask(0),
    fOutput(0)
{
/// Default constructor
  fLocalOutput[0] = fLocalOutput[1] = 0;

}
//----------------------------------------------------------------------
AliMUONRegionalTrigger::AliMUONRegionalTrigger(const AliMUONRegionalTrigger& theMUONRegionalTrig)
  : TObject(theMUONRegionalTrig),
    fId(theMUONRegionalTrig.fId),
    fLocalMask(theMUONRegionalTrig.fLocalMask),  
    fOutput(theMUONRegionalTrig.fOutput)           
{
/// Copy constructor (useful for TClonesArray)

  fLocalOutput[0] = theMUONRegionalTrig.fLocalOutput[0];
  fLocalOutput[1] = theMUONRegionalTrig.fLocalOutput[1];
 
}
//----------------------------------------------------------------------
AliMUONRegionalTrigger& AliMUONRegionalTrigger::operator=(const AliMUONRegionalTrigger& theMUONRegionalTrig)
{
/// Assigment operator;
/// equal operator (useful for non-pointer member in TClonesArray)

  if (this == &theMUONRegionalTrig)
    return *this;

  // base class assignement
  TObject::operator=(theMUONRegionalTrig);

  fId             = theMUONRegionalTrig.fId;
  fLocalMask      = theMUONRegionalTrig.fLocalMask;   
  fLocalOutput[0] = theMUONRegionalTrig.fLocalOutput[0];
  fLocalOutput[1] = theMUONRegionalTrig.fLocalOutput[1];
  fOutput         = theMUONRegionalTrig.fOutput;           

  return *this;
}



//----------------------------------------------------------------------
void AliMUONRegionalTrigger::Print(Option_t* opt) const
{
  //
  // Printing Regional Trigger information
  //
  TString sopt(opt);
  sopt.ToUpper();
 
  if ( sopt.Contains("FULL") ) { 

      printf("<AliMUONRegionalTrigger> Id %d localMask %d localOutputs %d %d output %d\n",
	     fId, fLocalMask, fLocalOutput[0], fLocalOutput[1], fOutput);

  }
}

