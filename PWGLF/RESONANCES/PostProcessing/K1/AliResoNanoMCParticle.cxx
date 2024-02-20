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

/* AliResoNanoMCParticle.h
 * A simple class to store MC particle information in a reduced format
 *
 *
 * Author: Bong-Hwi Lim
 *
 */

#include "AliResoNanoMCParticle.h"

ClassImp(AliResoNanoMCParticle);

AliResoNanoMCParticle::AliResoNanoMCParticle() : TNamed(),
                                                 fID(-1),
                                                 fEventID(-1),
                                                 fMotherID(-1),
                                                 fDaughter{-1, -1},
                                                 fPDGCode(-999),
                                                 fPt(-999),
                                                 fEta(-999),
                                                 fPhi(-999),
                                                 fRap(-999),
                                                 fCharge(-999),
                                                 fStatus(-999)
{
}

AliResoNanoMCParticle::~AliResoNanoMCParticle()
{
}
