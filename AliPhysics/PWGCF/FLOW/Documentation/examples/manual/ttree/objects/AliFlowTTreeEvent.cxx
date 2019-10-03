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

/* example event class which can be stored in a ttree 
 * see marcos in PWGCF/FLOW/Documentation/examples/manual/ttree/macros for
 * further explanation
 * author: redmer alexander bertens (rbertens@cern.ch)
 */


#include "AliFlowTTreeEvent.h"

ClassImp(AliFlowTTreeEvent)
AliFlowTTreeEvent::AliFlowTTreeEvent(): TObject(),
  fRun(-1),            
  fV0M(-1),
  fTRK(-1), 
  fZvtx(+999)
{
  // default constructor
}

//_________________________________________________________
AliFlowTTreeEvent::~AliFlowTTreeEvent(){}
