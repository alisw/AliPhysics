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

#include "AliMUONLocalStruct.h"

/// \class AliMUONLocalStruct
/// Local structure for trigger raw data.
/// The structure includes the information
///  about the x,y position of the 4 detection planes,
/// the trigger word (address, local decision, y trigger, y position, x deviation,
/// x position)
///
/// \author Christian Finck

/// \cond CLASSIMP
ClassImp(AliMUONLocalStruct)
/// \endcond

 const Int_t  AliMUONLocalStruct::fgkLength = 5;
 const Int_t  AliMUONLocalStruct::fgkScalerLength = 45;
 const UInt_t AliMUONLocalStruct::fgkEndOfLocal   = 0xCAFEFADE;
 const UInt_t AliMUONLocalStruct::fgkDisableWord  = 0x010CDEAD;
//___________________________________________
AliMUONLocalStruct::AliMUONLocalStruct()
  :  TObject(),
     fL0(0),   
     fHold(0), 
     fClk(0),   
     fLPtNTrig(0), 
     fHPtNTrig(0), 
     fLPtRTrig(0), 
     fHPtRTrig(0), 
     fLPtLTrig(0), 
     fHPtLTrig(0), 
     fLPtSTrig(0), 
     fHPtSTrig(0), 
     fEOS(0),         
     fReset(0)       
{
  ///
  /// ctor
  ///
  for (Int_t i = 0; i < 5; i++)
    fData[i] = 0;

  for (Int_t i = 0; i < 8*4; i++)
    fScaler[i] = 0;


}

//___________________________________________
AliMUONLocalStruct::AliMUONLocalStruct(const AliMUONLocalStruct& event)
  :  TObject(event),
     fL0(event.fL0),
     fHold(event.fHold),
     fClk(event.fClk),
     fLPtNTrig(event.fLPtNTrig),
     fHPtNTrig(event.fHPtNTrig),
     fLPtRTrig(event.fLPtRTrig),
     fHPtRTrig(event.fHPtRTrig),
     fLPtLTrig(event.fLPtLTrig),
     fHPtLTrig(event.fHPtLTrig),
     fLPtSTrig(event.fLPtSTrig),
     fHPtSTrig(event.fHPtSTrig),
     fEOS(event.fEOS),
     fReset(event.fReset)
{
  ///
  /// copy ctor
  ///
  for (Int_t i = 0; i < 5; i++)
    fData[i] = event.fData[i];

  for (Int_t i = 0; i < 8*4; i++)
    fScaler[i] = event.fScaler[i];


}
//___________________________________________
AliMUONLocalStruct& 
AliMUONLocalStruct::operator=(const AliMUONLocalStruct& event)
{
  /// 
  /// assignment operator
  ///

  if (this == &event) return *this;

  fL0       = event.fL0;
  fHold     = event.fHold;
  fClk      = event.fClk;
  fLPtNTrig = event.fLPtNTrig;
  fHPtNTrig = event.fHPtNTrig;
  fLPtRTrig = event.fLPtRTrig;
  fHPtRTrig = event.fHPtRTrig;
  fLPtLTrig = event.fLPtLTrig;
  fHPtLTrig = event.fHPtLTrig;
  fLPtSTrig = event.fLPtSTrig;
  fHPtSTrig = event.fHPtSTrig;
  fEOS      = event.fEOS;
  fReset    = event.fReset;

  for (Int_t i = 0; i < 5; i++)
    fData[i] = event.fData[i];

  for (Int_t i = 0; i < 8*4; i++)
    fScaler[i] = event.fScaler[i];

  return *this;
}

//___________________________________________
void AliMUONLocalStruct::SetScalersNumbers()
{
  /// set numbers for scaler events for local structure
  /// crasy numbers for scaler words, while no beam is coming
  ///

  fL0       = 1000;   
  fHold     = 100; 
  fClk      = 10000;  
  fLPtNTrig = 1; 
  fHPtNTrig = 1; 
  fLPtRTrig = 2; 
  fHPtRTrig = 2; 
  fLPtLTrig = 3; 
  fHPtLTrig = 3; 
  fLPtSTrig = 4; 
  fHPtSTrig = 4; 
  fEOS      = 0x2AA;         
  fReset    = 10;     

  for (Int_t i = 0; i < 8*4; i++)
    fScaler[i] = i;

}
