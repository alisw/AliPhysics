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

// Event Header for Hydjet generator
// Output generator parameters are accessable
// for the user through this interface.
// Author: Rafael Diaz Valdes
// (rafael.diaz.valdes@cern.ch)
//

#include "AliGenHydjetEventHeader.h"
ClassImp(AliGenHydjetEventHeader)

AliGenHydjetEventHeader::AliGenHydjetEventHeader():
    fNjet(0),
    fImpactParam(0),
    fNbcol(0),
    fNpart(0),
    fNpyt(0),
    fNhyd(0)
{
    // Constructor
}

AliGenHydjetEventHeader::AliGenHydjetEventHeader(const char* name):
    AliGenEventHeader(name),
    fNjet(0),
    fImpactParam(0),
    fNbcol(0),
    fNpart(0),
    fNpyt(0),
    fNhyd(0)
{
    // Copy Constructor
}
