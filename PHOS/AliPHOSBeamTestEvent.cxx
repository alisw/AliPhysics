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

//_________________________________________________________________________
// Class for PHOS Beam Test event header. 
// Contains PHOS specific information 
//  about triggers, cristal temperature etc.
//  This class attached as a branch to TreeE 
//  and used mainly to keep trigger pattern for current event.    
// 
//*-- Author : Maxim Volkov (RRC KI) & Dmitri Peressounko (RRC KI & SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSBeamTestEvent.h"

ClassImp(AliPHOSBeamTestEvent)


//____________________________________________________________________________ 
AliPHOSBeamTestEvent::AliPHOSBeamTestEvent():
  fBeamEnergy(0.),
  fPattern(0)
{
//Nothig should be set by default.
}

//____________________________________________________________________________ 
  AliPHOSBeamTestEvent::~AliPHOSBeamTestEvent()
{

}
