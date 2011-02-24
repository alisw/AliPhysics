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

//
// Class AliHFEdetPIDqa
// Base class for detector PID QA describing the interface to the PID QA
// manager, keeping also commom functionality. The following functions have
// to be implemented by the detector PID QA classes:
//   Initialize (basic initialization, i.e. histograms)
//   ProcessTrack (filling of the QA container)
// The base class provides the ESD/AOD PID object for all detector PID QA 
// classes
//
// Author:
//   Markus Fasel <M.Fasel@gsi.de>
//

#include "AliAODpidUtil.h"
#include "AliESDpid.h"

#include "AliHFEdetPIDqa.h"

ClassImp(AliHFEdetPIDqa)

//____________________________________________________________
AliHFEdetPIDqa::AliHFEdetPIDqa():
    TNamed()
  , fQAmanager(NULL)
{
  //
  // Dummy constructor
  //
}

//____________________________________________________________
AliHFEdetPIDqa::AliHFEdetPIDqa(const Char_t *name, const Char_t *title):
    TNamed(name, title)
  , fQAmanager(NULL)
{
  //
  // Default constructor
  //
}

//____________________________________________________________
AliHFEdetPIDqa::AliHFEdetPIDqa(const AliHFEdetPIDqa &o):
    TNamed(o)
  , fQAmanager(o.fQAmanager)
{
  //
  // Copy constructor
  //
}

//____________________________________________________________
AliHFEdetPIDqa &AliHFEdetPIDqa::operator=(const AliHFEdetPIDqa &o){
  //
  // Make assignment
  //
  TNamed::operator=(o);

  fQAmanager = o.fQAmanager;
  return *this;
}

