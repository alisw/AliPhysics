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

#include <TMath.h>

#include "AliMUONChamberTrigger.h"
#include "AliMUONResponseTrigger.h"
#include "AliMUONHit.h"
#include "AliMUON.h"
#include "AliMUONGeometryTransformer.h"
#include "AliLog.h"

//-----------------------------------------------------------------------------
/// \class AliMUONChamberTrigger
///
/// Implementation of AliMUONChamber for the trigger
///
/// \deprecated This class is to be deprecated.
///
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONChamberTrigger)
/// \endcond

//-------------------------------------------

AliMUONChamberTrigger::AliMUONChamberTrigger()
  : AliMUONChamber(),
    fkGeomTransformer(0)
{
/// Default constructor
}

//-------------------------------------------

AliMUONChamberTrigger:: ~AliMUONChamberTrigger()
{
/// Destructor
}

//-------------------------------------------

AliMUONChamberTrigger::AliMUONChamberTrigger(Int_t id,
                              const AliMUONGeometryTransformer* kGeometryTransformer) 
  : AliMUONChamber(id),
    fkGeomTransformer(kGeometryTransformer)
{
/// Constructor using chamber id
}






