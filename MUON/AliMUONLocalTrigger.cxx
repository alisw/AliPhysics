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
#include <assert.h>
#include "AliLog.h"
#include "AliMUONLocalStruct.h"

/// \class AliMUONLocalTrigger
/// Local Trigger algorithm data outputs
/// (contains local trigger decision and bit patterns)                \n
/// Add SetLocalStruct method for rawdata  (Ch. Finck)
/// \author Ph. Crochet

/// \cond CLASSIMP
ClassImp(AliMUONLocalTrigger)
/// \endcond

//----------------------------------------------------------------------
AliMUONLocalTrigger::AliMUONLocalTrigger()
  : TObject(), 
    fLoCircuit(0),
    fLoStripX(0),
    fLoDev(0),
    fLoStripY(0),
    fLoLpt(0),
    fLoHpt(0),
    
    fX1Pattern(0),
    fX2Pattern(0),
    fX3Pattern(0),
    fX4Pattern(0),
    
    fY1Pattern(0),
    fY2Pattern(0),
    fY3Pattern(0),
    fY4Pattern(0),
    
    fLoDecision(0),
    fDigits(0)
{
/// Default constructor
}
//----------------------------------------------------------------------
AliMUONLocalTrigger::AliMUONLocalTrigger(const AliMUONLocalTrigger& theMUONLocalTrig)
    : TObject(theMUONLocalTrig),
      fLoCircuit(theMUONLocalTrig.fLoCircuit),
      fLoStripX(theMUONLocalTrig.fLoStripX),
      fLoDev(theMUONLocalTrig.fLoDev),
      fLoStripY(theMUONLocalTrig.fLoStripY),
      fLoLpt(theMUONLocalTrig.fLoLpt),
      fLoHpt(theMUONLocalTrig.fLoHpt),
      
      fX1Pattern(theMUONLocalTrig.fX1Pattern),
      fX2Pattern(theMUONLocalTrig.fX2Pattern),
      fX3Pattern(theMUONLocalTrig.fX3Pattern),
      fX4Pattern(theMUONLocalTrig.fX4Pattern),
      
      fY1Pattern(theMUONLocalTrig.fY1Pattern),
      fY2Pattern(theMUONLocalTrig.fY2Pattern),
      fY3Pattern(theMUONLocalTrig.fY3Pattern),
      fY4Pattern(theMUONLocalTrig.fY4Pattern),
      
      fLoDecision(theMUONLocalTrig.fLoDecision),
      fDigits(theMUONLocalTrig.fDigits)
{
/// Copy constructor (useful for TClonesArray)

}
//----------------------------------------------------------------------
AliMUONLocalTrigger& AliMUONLocalTrigger::operator=(const AliMUONLocalTrigger& theMUONLocalTrig)
{
/// Assigment operator;
/// equal operator (useful for non-pointer member in TClonesArray)

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

  fX1Pattern  = theMUONLocalTrig.fX1Pattern;
  fX2Pattern  = theMUONLocalTrig.fX2Pattern;
  fX3Pattern  = theMUONLocalTrig.fX3Pattern;
  fX4Pattern  = theMUONLocalTrig.fX4Pattern;

  fY1Pattern  = theMUONLocalTrig.fY1Pattern;
  fY2Pattern  = theMUONLocalTrig.fY2Pattern;
  fY3Pattern  = theMUONLocalTrig.fY3Pattern;
  fY4Pattern  = theMUONLocalTrig.fY4Pattern;

  fLoDecision =  theMUONLocalTrig.fLoDecision;

  fDigits = theMUONLocalTrig.fDigits;

  return *this;
}

//----------------------------------------------------------------------
AliMUONLocalTrigger::AliMUONLocalTrigger(const Int_t* localtr, const TArrayI& digits) :
    fLoCircuit(localtr[0]),
    fLoStripX(localtr[1]),      
    fLoDev(localtr[2]),      
    fLoStripY(localtr[3]),          
    fLoLpt(localtr[4]),
    fLoHpt(localtr[5]),
     // keep on with this way
    fX1Pattern((UShort_t)localtr[7]),
    fX2Pattern((UShort_t)localtr[8]),
    fX3Pattern((UShort_t)localtr[9]),
    fX4Pattern((UShort_t)localtr[10]),
    
    fY1Pattern((UShort_t)localtr[11]),
    fY2Pattern((UShort_t)localtr[12]),
    fY3Pattern((UShort_t)localtr[13]),
    fY4Pattern((UShort_t)localtr[14]),
    
    fLoDecision(0),
    fDigits(digits)
{
/// Add a local trigger object
}

//----------------------------------------------------------------------
Char_t AliMUONLocalTrigger::GetLoDecision()
{
/// Get local decision 
/// from H(L)pt;
/// returns local trigger decision

  fLoDecision  = (fLoLpt & 0x3);
  fLoDecision |= (fLoHpt << 2) & 0xC;

  return fLoDecision;
}

//----------------------------------------------------------------------
void AliMUONLocalTrigger::GetDigit(
		Int_t i, Int_t& chamber, Int_t& cathode, Int_t& digit) const
{
/// Returns the i'th digit that fired this circuit.
/// The number of digits can be found with NumberOfDigits(), that is 
/// i is valid in the range [ 0 , NumberOfDigits() - 1 ]

	Int_t digitnumber = fDigits[i];
	DecodeDigitNumber(digitnumber, chamber, cathode, digit);
}

//----------------------------------------------------------------------
Int_t AliMUONLocalTrigger::EncodeDigitNumber(Int_t chamber, Int_t cathode, 
					     Int_t digit)
{
/// Encodes a 32-bit digit number from digit information to be stored
/// in internal integer arrays. Note that the value of the digit parameter
/// can not be larger than 0x07FFFFFF.

	assert( 0 <= cathode && cathode <= 1 );
	assert( 0 <= chamber && chamber <= 13 );

	if ((digit & 0xF8000000) != 0)
	{
		AliErrorGeneral("AliMUONLocalTrigger", Form(
			"Digit index value is to large: 0x%8.8X. Maximum supported value is 0x07FFFFFF.",
			cathode, chamber, digit
		));
		return -1;
	};

	return ((cathode & 0x1) << 31) | ((chamber & 0xF) << 27) | digit;
}

//----------------------------------------------------------------------
void AliMUONLocalTrigger::DecodeDigitNumber(Int_t digitnumber, Int_t& chamber, 
					    Int_t& cathode, Int_t& digit)
{
/// Decodes a digit number into information about the digit.
/// One can subsequently fetch the digit with
/// AliMUONDataInterface::Digit(chamber, cathode, digit)

	cathode = (digitnumber >> 31) & 0x1;
	chamber = (digitnumber >> 27) & 0xF;
	digit = digitnumber & 0x7FFFFFF;
}

//----------------------------------------------------------------------
void AliMUONLocalTrigger::SetLocalStruct(Int_t loCircuit, AliMUONLocalStruct& localStruct)
{
/// Set local trigger info from rawdata localStruct

  // set id'
  SetLoCircuit(loCircuit);

  // set X, Y dev  
  SetLoStripX((Int_t)localStruct.GetXPos());
  SetLoStripY((Int_t)localStruct.GetYPos());
  SetLoDev((Int_t)localStruct.GetXDev());
 
  // set L(H)pt
  SetLoLpt(localStruct.GetLpt());
  SetLoHpt(localStruct.GetHpt());

  // set pattern X
  SetX1Pattern(localStruct.GetX1());
  SetX2Pattern(localStruct.GetX2());
  SetX3Pattern(localStruct.GetX3());
  SetX4Pattern(localStruct.GetX4());

  // set pattern Y
  SetY1Pattern(localStruct.GetY1());
  SetY2Pattern(localStruct.GetY2());
  SetY3Pattern(localStruct.GetY3());
  SetY4Pattern(localStruct.GetY4());

}

//----------------------------------------------------------------------
void AliMUONLocalTrigger::Print(Option_t* opt) const
{
  //
  // Printing Local Trigger information
  //
  TString sopt(opt);
  sopt.ToUpper();
 
  if ( sopt.Contains("FULL") ) { 

      printf("<AliMUONLocalTrigger> Circuit %d StripX %d Dev %d StripY %d Lpt %d Hpt %d \n",LoCircuit(),LoStripX(),LoDev(),LoStripY(),LoLpt(),LoHpt());

  }
}

