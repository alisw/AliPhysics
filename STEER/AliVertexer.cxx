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

#include "AliLog.h"
#include <AliVertexer.h>

ClassImp(AliVertexer)

//////////////////////////////////////////////////////////////////////
// Base class for primary vertex reconstruction                     //
// AliESDVertexer is a class for full 3D primary vertex finding     //
// derived classes: AliITSVertexer                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliVertexer::AliVertexer() {
  // Default Constructor

    fCurrentVertex  = 0;
    SetFirstEvent(0);
    SetLastEvent(0);
}


//______________________________________________________________________
AliVertexer::AliVertexer(const AliVertexer &vtxr) : TObject(vtxr) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  AliFatal("Copy constructor not allowed");
}

//______________________________________________________________________
AliVertexer& AliVertexer::operator=(const AliVertexer& /* vtxr */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  AliFatal("Assignment operator not allowed");
  return *this;
}

//______________________________________________________________________
AliVertexer::~AliVertexer() {
  // Default Destructor
  // The objects pointed by the following pointers are not owned
  // by this class and are not deleted

    fCurrentVertex  = 0;
}

//______________________________________________________________________
void AliVertexer::SetDebug(Int_t debug)
{
  AliWarning("Don't use this method any more, use AliDebug instead");
  fDebug = debug;
}
