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

#include "AliMUONLocalTrigger.h"

ClassImp(AliMUONLocalTrigger);
//----------------------------------------------------------------------
AliMUONLocalTrigger::AliMUONLocalTrigger()
  : TObject()
{
// constructor
  fLoCircuit = 0;
  fLoStripX  = 0;
  fLoDev     = 0;
  fLoStripY  = 0;
  fLoLpt     = 0;
  fLoHpt     = 0;
  fLoApt     = 0;
}
//----------------------------------------------------------------------
AliMUONLocalTrigger::AliMUONLocalTrigger(const AliMUONLocalTrigger& theMUONLocalTrig)
  : TObject(theMUONLocalTrig)
{
// copy constructor (useful for TClonesArray)
  fLoCircuit = theMUONLocalTrig.fLoCircuit;
  fLoStripX  = theMUONLocalTrig.fLoStripX;         
  fLoDev     = theMUONLocalTrig.fLoDev;           
  fLoStripY  = theMUONLocalTrig.fLoStripY;           
  fLoLpt     = theMUONLocalTrig.fLoLpt;
  fLoHpt     = theMUONLocalTrig.fLoHpt;
  fLoApt     = theMUONLocalTrig.fLoApt;
}
//----------------------------------------------------------------------
AliMUONLocalTrigger& AliMUONLocalTrigger::operator=(const AliMUONLocalTrigger& theMUONLocalTrig)
{
// equal operator (useful for non-pointer member in TClonesArray)
  if (this == &theMUONLocalTrig)
    return *this;

  // base class assignement
  TObject::operator=(theMUONLocalTrig);

  fLoCircuit = theMUONLocalTrig.fLoCircuit;
  fLoStripX  = theMUONLocalTrig.fLoStripX;         
  fLoDev     = theMUONLocalTrig.fLoDev;           
  fLoStripY  = theMUONLocalTrig.fLoStripY;           
  fLoLpt     = theMUONLocalTrig.fLoLpt;
  fLoHpt     = theMUONLocalTrig.fLoHpt;
  fLoApt     = theMUONLocalTrig.fLoApt;

  return *this;
}

//----------------------------------------------------------------------
AliMUONLocalTrigger::AliMUONLocalTrigger(Int_t *localtr)
{
// add a local trigger object 
  fLoCircuit = localtr[0];
  fLoStripX  = localtr[1];         
  fLoDev     = localtr[2];           
  fLoStripY  = localtr[3];           
  fLoLpt     = localtr[4];
  fLoHpt     = localtr[5];
  fLoApt     = localtr[6];
}

