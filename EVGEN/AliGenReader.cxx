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
//
// Interface for reading events from files.
// Realisations of this interface have to be used with AliGenExFile.
// NextEvent() loops over events 
// and NextParticle() loops over particles. 
// Author: andreas.morsch@cern.ch

#include "AliGenReader.h"
ClassImp(AliGenReader)


AliGenReader& AliGenReader::operator=(const  AliGenReader& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

void AliGenReader::Copy(AliGenReader&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}



