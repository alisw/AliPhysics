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

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides the minimum-size TOF hit info       //
//                                                           //
//   author: Roberto Preghenella (R+)                        //
//           preghenella@bo.infn.it                          //
//                                                           //
///////////////////////////////////////////////////////////////

#include "AliTOFHitField.h"

ClassImp(AliTOFHitField)

//-----------------------------------------------------------------------------

AliTOFHitField::AliTOFHitField() : 
  fIndex(0),
  fTimeBin(0),
  fTOTBin(0),
  fDeltaBC(0),
  fL0L1Latency(0)
{
  /*
   * default constructor
   */
}

//-----------------------------------------------------------------------------

AliTOFHitField::~AliTOFHitField()
{
  /*
   * default destructor
   */
}

//-----------------------------------------------------------------------------

AliTOFHitField::AliTOFHitField(const AliTOFHitField &source) :
  fIndex(source.fIndex),
  fTimeBin(source.fTimeBin),
  fTOTBin(source.fTOTBin),
  fDeltaBC(source.fDeltaBC),
  fL0L1Latency(source.fL0L1Latency)
{ 
  /*
   * default constructor
   */
}

//-----------------------------------------------------------------------------

AliTOFHitField & 
AliTOFHitField::operator=(const AliTOFHitField &source) 
{ 
  /*
   * operator=
   */

  if(this==&source) return *this;
  fIndex = source.fIndex;
  fTimeBin = source.fTimeBin;
  fTOTBin = source.fTOTBin;
  fDeltaBC = source.fDeltaBC;
  fL0L1Latency =source.fL0L1Latency;
  return *this;
}

