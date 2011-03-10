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


/////////////////////////////////////////////////////////////////////////////////
// This class sets the local coordinates via a specific setter. Needed because //
// the AliGeomManager class can not be used for the upgrade code at this stage //
/////////////////////////////////////////////////////////////////////////////////

#include <AliITSRecPointU.h>
ClassImp(AliITSRecPointU)
//_____________________________________________________________
AliITSRecPointU::AliITSRecPointU():
    AliITSRecPoint(),
    fModule(0)
{
 //
 // Default constructor
 // 
}
//_____________________________________________________________
AliITSRecPointU::AliITSRecPointU(const AliITSRecPointU& pt):
    AliITSRecPoint(pt),
    fModule(pt.fModule)
{
  //
  // Copy constructor
  //
}
//______________________________________________________________________
AliITSRecPointU& AliITSRecPointU::operator=(const AliITSRecPointU& source)
{
  //
  // Assignment operator (as in AliITSRecPoint)
  //

  this->~AliITSRecPointU();
  new(this) AliITSRecPointU(source);
  return *this;

}
