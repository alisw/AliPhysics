/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//     Event cuts base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include "AliVCuts.h"

ClassImp(AliVCuts)

//______________________________________________________________________________
AliVCuts::AliVCuts() : 
  TNamed("Cuts","") { } // default constructor 

AliVCuts::AliVCuts(const char* name, const char* title) :
  TNamed(name, title) { } 
//______________________________________________________________________________
AliVCuts::AliVCuts(const AliVCuts& cuts) :
  TNamed(cuts) {} // Copy constructor

//______________________________________________________________________________
AliVCuts& AliVCuts::operator=(const AliVCuts& cuts)
{
  // Assignment operator
  if(this!=&cuts) { 
    TNamed::operator=(cuts);
  }
  return *this;
}
