/**************************************************************************
* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

#include <TTree.h>
#include <TFile.h>
#include "AliLog.h"
#include "AliCSTrackCutsBase.h"

/// \file AliCSTrackCutsBase.cxx
/// \brief Implementation of abstract base track cuts class within the correlation studies analysis

/// Default constructor for serialization
AliCSTrackCutsBase::AliCSTrackCutsBase() :
    AliCSAnalysisCutsBase()
{
}

/// Constructor
/// Allocates the needed memory for the number of cuts to support
/// \param nCuts the number of cuts to support
/// \param nParams the number of cuts parameters
/// \param name name of the event cuts
/// \param title title of the event cuts
AliCSTrackCutsBase::AliCSTrackCutsBase(Int_t nCuts, Int_t nParams, const char *name, const char *title) :
    AliCSAnalysisCutsBase(nCuts, nParams, name, title)
{
}

/// Destructor
AliCSTrackCutsBase::~AliCSTrackCutsBase()
{
}


/// \cond CLASSIMP
ClassImp(AliCSTrackCutsBase);
/// \endcond
