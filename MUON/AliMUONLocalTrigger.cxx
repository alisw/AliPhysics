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

ClassImp(AliMUONLocalTrigger)
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

  fX1Pattern  = 0;
  fX2Pattern  = 0;
  fX3Pattern  = 0;
  fX4Pattern  = 0;

  fY1Pattern  = 0;
  fY2Pattern  = 0;
  fY3Pattern  = 0;
  fY4Pattern  = 0;

  fLoDecision = 0;

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

  fX1Pattern  = theMUONLocalTrig.fX1Pattern;
  fX2Pattern  = theMUONLocalTrig.fX2Pattern;
  fX3Pattern  = theMUONLocalTrig.fX3Pattern;
  fX4Pattern  = theMUONLocalTrig.fX4Pattern;

  fY1Pattern  = theMUONLocalTrig.fY1Pattern;
  fY2Pattern  = theMUONLocalTrig.fY2Pattern;
  fY3Pattern  = theMUONLocalTrig.fY3Pattern;
  fY4Pattern  = theMUONLocalTrig.fY4Pattern;

  fLoDecision =  theMUONLocalTrig.fLoDecision;
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

  fX1Pattern  = theMUONLocalTrig.fX1Pattern;
  fX2Pattern  = theMUONLocalTrig.fX2Pattern;
  fX3Pattern  = theMUONLocalTrig.fX3Pattern;
  fX4Pattern  = theMUONLocalTrig.fX4Pattern;

  fY1Pattern  = theMUONLocalTrig.fY1Pattern;
  fY2Pattern  = theMUONLocalTrig.fY2Pattern;
  fY3Pattern  = theMUONLocalTrig.fY3Pattern;
  fY4Pattern  = theMUONLocalTrig.fY4Pattern;

  fLoDecision =  theMUONLocalTrig.fLoDecision;

  return *this;
}

//----------------------------------------------------------------------
AliMUONLocalTrigger::AliMUONLocalTrigger(Int_t* localtr)
{
// add a local trigger object 
  fLoCircuit = localtr[0];
  fLoStripX  = localtr[1];         
  fLoDev     = localtr[2];           
  fLoStripY  = localtr[3];           
  fLoLpt     = localtr[4];
  fLoHpt     = localtr[5];
  fLoApt     = localtr[6];

  // keep on with this way
  fX1Pattern = (UShort_t)localtr[7];
  fX2Pattern = (UShort_t)localtr[8];
  fX3Pattern = (UShort_t)localtr[9];
  fX4Pattern = (UShort_t)localtr[10];

  fY1Pattern = (UShort_t)localtr[11];
  fY2Pattern = (UShort_t)localtr[12];
  fY3Pattern = (UShort_t)localtr[13];
  fY4Pattern = (UShort_t)localtr[14];
}
//----------------------------------------------------------------------
Char_t AliMUONLocalTrigger::GetLoDecision()
{
  fLoDecision  = (fLoLpt & 0x3);
  fLoDecision |= (fLoHpt << 2) & 0xC;

  return fLoDecision;
}
