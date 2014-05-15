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

/* $Id: AliGenEventHeaderTunedPbPb.cxx 50615 2013-08-19 23:19:19Z fnoferin $ */

//---------------------------------------------------------------------
// Event header base class for generator. 
// Stores as a minimum the date, run number, event number,
// number of particles produced  
// and the impact parameter.
// + Psi_2 Psi_3 Psi4
// Author: fnoferin@cern.ch
//---------------------------------------------------------------------

#include "AliGenEventHeaderTunedPbPb.h"
ClassImp(AliGenEventHeaderTunedPbPb)


//_______________________________________________________________________
AliGenEventHeaderTunedPbPb::AliGenEventHeaderTunedPbPb():
  AliGenEventHeader(),
  fCentrality(0.),
  fPsi2(0.),
  fPsi3(0.),
  fPsi4(0.)
{
  //
  // Constructor
  //
}

//_______________________________________________________________________
AliGenEventHeaderTunedPbPb::AliGenEventHeaderTunedPbPb(const char * name):
  AliGenEventHeader(name),
  fCentrality(0.),
  fPsi2(0.),
  fPsi3(0.),
  fPsi4(0.)
{
  //
  // Constructor
  //
}

