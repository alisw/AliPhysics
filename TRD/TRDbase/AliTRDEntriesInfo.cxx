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

/* $Id: AliTRDEntriesInfo.cxx 27946 2008-08-13 15:26:24Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Calibration base class for a single ROC                                  //
//  Contains one UShort_t value per pad                                      //
//  However, values are set and get as float, there are stored internally as //
//  (UShort_t) value * 10000                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDEntriesInfo.h"

ClassImp(AliTRDEntriesInfo)

//_____________________________________________________________________________
AliTRDEntriesInfo::AliTRDEntriesInfo()
  :AliTRDUshortInfo()
{
  //
  // Default constructor
  //

}
//_____________________________________________________________________________
AliTRDEntriesInfo::AliTRDEntriesInfo(Int_t n)
  :AliTRDUshortInfo(n)
{
  //
  // Constructor that initializes a given size
  //
 
}
//_____________________________________________________________________________
AliTRDEntriesInfo::AliTRDEntriesInfo(const AliTRDEntriesInfo &c)
  :AliTRDUshortInfo(c)
{
  //
  // AliTRDEntriesInfo copy constructor
  //
  
}
//_____________________________________________________________________________
AliTRDEntriesInfo::~AliTRDEntriesInfo()
{
  //
  // AliTRDEntriesInfo destructor
  //

  
}
//_____________________________________________________________________________
AliTRDEntriesInfo &AliTRDEntriesInfo::operator=(const AliTRDEntriesInfo &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDEntriesInfo &) c).Copy(*this);
  return *this;

}
//___________________________________________________________________________________
Int_t AliTRDEntriesInfo::GetSum() const
{
  //
  // Calculate the sum of entries
  //

  Int_t total = 0;
  
  for(Int_t k = 0; k < fSize; k++){
    total += fData[k];
  }


  return total;
  
}
//____________________________________________________________________________________________
Bool_t AliTRDEntriesInfo::TestAdd(const AliTRDEntriesInfo * info)
{
  //
  // add values 
  //
  for (Int_t  idata = 0; idata< fSize; idata++){
    if((At(idata)+info->At(idata)) > 65535) return kFALSE;
  }
  return kTRUE;
}
//____________________________________________________________________________________________
void AliTRDEntriesInfo::Add(const AliTRDEntriesInfo * info)
{
  //
  // add values 
  //
  for (Int_t  idata = 0; idata< fSize; idata++){
    fData[idata] += info->At(idata);    
  }
}
//____________________________________________________________________________________________
void AliTRDEntriesInfo::AddIf(const AliTRDEntriesInfo * info)
{
  //
  // add values 
  //
  for (Int_t  idata = 0; idata< fSize; idata++){
    if(((fData[idata]+info->At(idata)) <= 65535) && ((fData[idata]+info->At(idata)) >= 0)) fData[idata] += info->At(idata);    
  }
}
