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
 
#include "AliTender.h"
#include "AliTenderSupply.h"

ClassImp(AliTenderSupply)

//______________________________________________________________________________
AliTenderSupply::AliTenderSupply()
                :TNamed(),
                 fTender(NULL)
{
// Dummy constructor
}

//______________________________________________________________________________
AliTenderSupply::AliTenderSupply(const char* name, const AliTender *tender)
                :TNamed(name, "ESD analysis tender car"),
                 fTender(tender)
{
// Default constructor
}

//______________________________________________________________________________
AliTenderSupply::AliTenderSupply(const AliTenderSupply &other)
                :TNamed(other),
                 fTender(other.fTender)
                 
{
// Copy constructor
}

//______________________________________________________________________________
AliTenderSupply::~AliTenderSupply()
{
// Destructor
}

//______________________________________________________________________________
AliTenderSupply& AliTenderSupply::operator=(const AliTenderSupply &other)
{
// Assignment
   if (&other == this) return *this;
   TNamed::operator=(other);
   fTender = other.fTender;
   return *this;
}
