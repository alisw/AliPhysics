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

#include "AliGenReader.h"
ClassImp(AliGenReader)


AliGenReader& AliGenReader::operator=(const  AliGenReader& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

void AliGenReader::RewindEvent()
{
  // Go back to the first particle of the event.
  // Need to be implemented in the implementation classes. Interface dies.
  Fatal("AliGenReader::RewindEvent","\nMethod RewindEvent not present in the implementation class.\n");
}

void AliGenReader::Copy(AliGenReader&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}



