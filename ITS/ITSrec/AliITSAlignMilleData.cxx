/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------------
/// \class AliITSAlignMilleData
/// Alignment class for the ALICE ITS detector
///
///
/// \author M. Lunardon 
//-----------------------------------------------------------------------------

#include "AliITSAlignMilleData.h"

/// \cond CLASSIMP
ClassImp(AliITSAlignMilleData)
/// \endcond

AliITSAlignMilleData::AliITSAlignMilleData() 
: TObject(),
  fMeasX(0),
  fSigmaX(0), 
  fMeasZ(0),
  fSigmaZ(0) 
{
  // Default constructor
  for (Int_t i=0; i<ITSMILLENLOCAL; i++) {
    fIdxlocX[i]=0;
    fDerlocX[i]=0.0;
    fIdxlocZ[i]=0;
    fDerlocZ[i]=0.0;
  }
  for (Int_t i=0; i<ITSMILLENPARCH; i++) {
    fIdxgloX[i]=0;
    fDergloX[i]=0.0;
    fIdxgloZ[i]=0;
    fDergloZ[i]=0.0;
  }
}

AliITSAlignMilleData::~AliITSAlignMilleData() {
  /// Destructor

}

