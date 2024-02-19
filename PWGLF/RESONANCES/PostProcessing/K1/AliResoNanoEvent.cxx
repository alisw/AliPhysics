/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/* AliResoNanoEvent.h
 * Header file for the AliResoNanoEvent class, which is a simple class
 * to store event information in a reduced format
 *
 * Author: Bong-Hwi Lim
 *
 */
#include "AliResoNanoEvent.h"
class AliResoNanoEvent;

AliResoNanoEvent::AliResoNanoEvent() : fEventID(-1),
                                       fCentrality(-1),
                                       fVtxX(-999),
                                       fVtxY(-999),
                                       fVtxZ(-999),
                                       fVtxZErr(-999),
                                       fVtxMCX(-999),
                                       fVtxMCY(-999),
                                       fVtxMCZ(-999)
{
}

AliResoNanoEvent::~AliResoNanoEvent()
{
}
