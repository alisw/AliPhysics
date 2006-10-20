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

// $Id$

#include "AliMUONVSubprocessor.h"

/// \class AliMUONVSubprocessor
///
/// Base class for a shuttle sub-task for MUON (either TRK or TRG)
///
/// It allows to organize more cleanly the various calibration tasks
/// to be performed within the MUON shuttle's preprocessors
///
/// \author Laurent Aphecetche
///

/// \cond CLASSIMP
ClassImp(AliMUONVSubprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONVSubprocessor::AliMUONVSubprocessor(AliMUONPreprocessor* master,
                                           const char* name,
                                           const char* title)
: TNamed(name,title), fMaster(master)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONVSubprocessor::~AliMUONVSubprocessor()
{
  /// dtor
}

//_____________________________________________________________________________
void
AliMUONVSubprocessor::Initialize(Int_t /*run*/, 
                                 UInt_t /*startTime*/, 
                                 UInt_t /*endTime*/)
{
  /// optional
}

