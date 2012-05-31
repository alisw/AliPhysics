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

/* $Id: AliTRDPhInfo.cxx 27946 2008-08-13 15:26:24Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Calibration base class for a single ROC                                  //
//  Contains one UShort_t value per pad                                      //
//  However, values are set and get as float, there are stored internally as //
//  (UShort_t) value * 10000                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDPhInfo.h"

ClassImp(AliTRDPhInfo)

//_____________________________________________________________________________
AliTRDPhInfo::AliTRDPhInfo()
  :AliTRDUshortInfo()
{
  //
  // Default constructor
  //

}
//_____________________________________________________________________________
AliTRDPhInfo::AliTRDPhInfo(Int_t n)
  :AliTRDUshortInfo(n)
{
  //
  // Constructor that initializes a given size
  //
 
}
//_____________________________________________________________________________
AliTRDPhInfo::AliTRDPhInfo(const AliTRDPhInfo &c)
  :AliTRDUshortInfo(c)
{
  //
  // AliTRDPhInfo copy constructor
  //
  
}
//_____________________________________________________________________________
AliTRDPhInfo::~AliTRDPhInfo()
{
  //
  // AliTRDPhInfo destructor
  //

  
}
//_____________________________________________________________________________
AliTRDPhInfo &AliTRDPhInfo::operator=(const AliTRDPhInfo &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDPhInfo &) c).Copy(*this);
  return *this;

}

