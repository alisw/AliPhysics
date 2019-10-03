/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////////////
//
// Compatibility wrapper for AliCanvas
//
// Author: Jochen Klein <jochen.klein@cern.ch>
//
////////////////////////////////////////////////////////////////////////////

#include "TError.h"

#include "AliFigure.h"

AliFigure::AliFigure(const char* name, const char* title, Int_t ww, Int_t wh) :
  AliCanvas(name, title, ww, wh)
{
  ::Obsolete("AliFigure -> AliCanvas", "vAN-20150710", "vAN-20151231");
}
