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
#include "AliESDVertex.h"
#include "AliVertexer.h"
#include "AliMultiplicity.h"

ClassImp(AliVertexer)

//////////////////////////////////////////////////////////////////////
// Base class for primary vertex reconstruction                     //
// AliESDVertexer is a class for full 3D primary vertex finding     //
// derived classes: AliITSVertexer                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliVertexer::AliVertexer() :
  fCurrentVertex(NULL),
  fMult(NULL)
{
  //
  // Default Constructor
  //
  SetVtxStart(0.,0.,0.);
  for(Int_t i=0;i<6;i++)fNominalCov[i]=0.;
}

//______________________________________________________________________
AliVertexer::~AliVertexer() {
  // Default Destructor

  if(fMult) delete fMult;
}


//---------------------------------------------------------------------------
void  AliVertexer::SetVtxStart(AliESDVertex *vtx) 
{ 
//
// Set initial vertex knowledge
//
  vtx->GetXYZ(fNominalPos);
  vtx->GetCovMatrix(fNominalCov);
  return; 
}
