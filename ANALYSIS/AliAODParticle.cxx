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

/////////////////////////////////////////////////////////////
//
// base class for AOD particles
//
// Most of the methods defined in this base class (or even all of them)
// should be redefined in derived classes. Some methods like Pt() or
// Vx() have default implementation, but they are not optimised for
// performance. It is likely that they can be (and should be) implemented
// in the derived classes in a way which gives faster access.
//
// Algorithms and analysis classes should as much as possible use the
// interface of the base class instead of special methods of derived
// classes. This allows to run the code on all types of particles.
//
/////////////////////////////////////////////////////////////

#include "AliAODParticle.h"

Int_t AliAODParticle::fgDebug = 0;

ClassImp(AliAODParticle)

