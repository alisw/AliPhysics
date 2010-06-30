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
#include "AliLog.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONRawStreamTriggerHP.h"
#include <Riostream.h>
#include <TArrayS.h>

//-----------------------------------------------------------------------------
/// \class AliMUONLocalTrigger
/// Local Trigger algorithm data outputs
/// (contains local trigger decision and bit patterns)                \n
/// Add SetLocalStruct method for rawdata  (Ch. Finck)
/// \author Ph. Crochet
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONLocalTrigger)
/// \endcond

//----------------------------------------------------------------------
AliMUONLocalTrigger::AliMUONLocalTrigger()
  : TObject(), 
    fLoCircuit(0),
    fLoStripX(0),
    fLoDev(0),
    fLoSdev(1),
    fLoTrigY(1),
    fLoStripY(15),
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

    fHitPatternFromResponse(0xFF),
    fTriggerWithoutChamber(0)
{
/// Default constructor
}
//----------------------------------------------------------------------
AliMUONLocalTrigger::AliMUONLocalTrigger(const AliMUONLocalTrigger& theMUONLocalTrig)
    : TObject(theMUONLocalTrig),
      fLoCircuit(theMUONLocalTrig.fLoCircuit),
      fLoStripX(theMUONLocalTrig.fLoStripX),
      fLoDev(theMUONLocalTrig.fLoDev),
      fLoSdev(theMUONLocalTrig.fLoSdev),
      fLoTrigY(theMUONLocalTrig.fLoTrigY),
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

      fHitPatternFromResponse(theMUONLocalTrig.fHitPatternFromResponse),
      fTriggerWithoutChamber(theMUONLocalTrig.fTriggerWithoutChamber)
{
/// Copy constructor (useful for TClonesArray)

}

//----------------------------------------------------------------------
AliMUONLocalTrigger::~AliMUONLocalTrigger()
{
/// Destructor
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
  fLoSdev    = theMUONLocalTrig.fLoSdev;           
  fLoTrigY   = theMUONLocalTrig.fLoTrigY;           
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

  fHitPatternFromResponse = theMUONLocalTrig.fHitPatternFromResponse;
  fTriggerWithoutChamber = theMUONLocalTrig.fTriggerWithoutChamber;

  return *this;
}


//----------------------------------------------------------------------
Char_t AliMUONLocalTrigger::GetLoDecision() const
{
/// Get local decision 
/// from H(L)pt;
/// returns local trigger decision

  Char_t rv = (fLoLpt & 0x3);
  rv |= (fLoHpt << 2) & 0xC;

  return rv;
}

//___________________________________________
void AliMUONLocalTrigger::GetXPattern(TArrayS& array) const
{
    /// return array of X pattern
    Short_t vec[4] = {GetX1Pattern(), GetX2Pattern(), GetX3Pattern(), GetX4Pattern()};
    array.Set(4, vec);
}

//___________________________________________
void AliMUONLocalTrigger::GetYPattern(TArrayS& array) const
{
    /// return array of Y pattern
    Short_t vec[4] = {GetY1Pattern(), GetY2Pattern(), GetY3Pattern(), GetY4Pattern()};
    array.Set(4, vec);
}

//___________________________________________
Bool_t
AliMUONLocalTrigger::IsNull() const
{
  /// Whether or not this card has something usefull to say or not
  return ( fX1Pattern == 0 &&
           fX2Pattern == 0 &&
           fX3Pattern == 0 &&
           fX4Pattern == 0 &&
           fY1Pattern == 0 &&
           fY2Pattern == 0 &&
           fY3Pattern == 0 &&
           fY4Pattern == 0 );          
}

//----------------------------------------------------------------------
void AliMUONLocalTrigger::SetLocalStruct(Int_t loCircuit, AliMUONLocalStruct& localStruct)
{
/// Set local trigger info from rawdata localStruct

  // set id'
  SetLoCircuit(loCircuit);

  // set X, Y, dev, Sdev and TrigY
  SetLoStripX((Int_t)localStruct.GetXPos());
  SetLoStripY((Int_t)localStruct.GetYPos());
  SetLoDev((Int_t)localStruct.GetXDev());
  SetLoSdev((Int_t)localStruct.GetSXDev());
  SetLoTrigY((Int_t)localStruct.GetTrigY());
 
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
void AliMUONLocalTrigger::SetLocalStruct(Int_t loCircuit, const AliMUONRawStreamTriggerHP::AliLocalStruct& localStruct)
{
/// Set local trigger info from rawdata localStruct (new raw reader)

  // set id'
  SetLoCircuit(loCircuit);

  // set X, Y, dev, Sdev and TrigY
  SetLoStripX((Int_t)localStruct.GetXPos());
  SetLoStripY((Int_t)localStruct.GetYPos());
  SetLoDev((Int_t)localStruct.GetXDev());
  SetLoSdev((Int_t)localStruct.GetSXDev());
  SetLoTrigY((Int_t)localStruct.GetTrigY());
 
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

namespace
{
  const char* AsString(Int_t t)
  {
    switch (t)
    {
      case 0:
        return "no";
        break;
      case 1:
        return "minus";
        break;
      case 2:
        return "plus";
        break;
      case 3:
        return "undef";
        break;
      default:
        return "";
        break;
    }
  }
}

//----------------------------------------------------------------------
void AliMUONLocalTrigger::Print(Option_t* opt) const
{
/// Printing Local Trigger information

  TString sopt(opt);
  sopt.ToUpper();

  cout << Form("Circuit %3d Decision %2d StripX %2d Dev %2d(%1d) StripY %2d Lpt %6s Hpt %6s",
               LoCircuit(), GetLoDecision(),
               LoStripX(), LoDev(), LoSdev(), LoStripY(),
               AsString(LoLpt()),AsString(LoHpt()),IsNull()) << endl;
  
  if ( sopt.Contains("FULL") ) { 

    cout << Form("Xpatterns = 0x %04x %04x %04x %04x",
                 fX1Pattern,fX2Pattern,fX3Pattern,fX4Pattern) << endl;
    cout << Form("Ypatterns = 0x %04x %04x %04x %04x",
                 fY1Pattern,fY2Pattern,fY3Pattern,fY4Pattern) << endl;
  }
}

//----------------------------------------------------------------------
Int_t AliMUONLocalTrigger::GetDeviation() const
{
  /// return deviation
  
  Int_t deviation = LoDev(); 
  Int_t sign = 0;
  if ( !LoSdev() &&  deviation ) sign=-1;
  if ( !LoSdev() && !deviation ) sign= 0;
  if (  LoSdev() == 1 )          sign=+1;
  deviation *= sign;
  deviation += 15;
  return deviation;
}

//----------------------------------------------------------------------
void AliMUONLocalTrigger::SetDeviation(Int_t deviation)
{
  /// set LoDev and LoSDev according to deviation
  
  deviation -= 15;
  if (deviation > 0) {
    SetLoDev(deviation);
    SetLoSdev(1);
  } else {
    SetLoDev(-deviation);
    SetLoSdev(0);
  }
}

//----------------------------------------------------------------------
const char*
AliMUONLocalTrigger::GetName() const
{
/// Generate name

  return Form("LocalBoard%3d",LoCircuit());
}


//----------------------------------------------------------------------
Bool_t AliMUONLocalTrigger::IsTrigX() const
{
/// Trigger response X strips
  Bool_t xTrig;
  if ( LoSdev()==1 && LoDev()==0 && 
       LoStripX()==0) xTrig=kFALSE; // no trigger in X
  else xTrig = kTRUE;                       // trigger in X
  return xTrig;
}


//----------------------------------------------------------------------
Bool_t AliMUONLocalTrigger::IsTrigY() const
{
/// Trigger response Y strips
  Bool_t yTrig;
  if ( LoTrigY()==1 && 
       LoStripY()==15 ) yTrig = kFALSE; // no trigger in Y
  else yTrig = kTRUE;                          // trigger in Y
  return yTrig;
}
