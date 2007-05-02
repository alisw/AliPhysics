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

// -------------------------
// Class AliMUONDataLoader
// -------------------------
// Implements AliLoader for MUON subsystem
// Author: Gines Martinez

#include "AliMUONSimLoader.h"
#include "AliMUONConstants.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONSimLoader)
/// \endcond

//___________________________________________________________________
AliMUONSimLoader::AliMUONSimLoader()
  : AliLoader(),
    fMUONData(0)
{
/// Default constructor
}
//_______________________________________________________________________________
AliMUONSimLoader::AliMUONSimLoader(const Char_t* detname,const Char_t* eventfoldername)
  : AliLoader(detname,eventfoldername),
    fMUONData(0)
{
/// Standard contructor with name of the top folder of the tree
}
//_______________________________________________________________________________
AliMUONSimLoader::AliMUONSimLoader(const Char_t * detname,TFolder* eventfolder)
  : AliLoader(detname,eventfolder),
    fMUONData(0)
{
/// Constructor
}

//_______________________________________________________________________________
AliMUONSimLoader::~AliMUONSimLoader()
{
/// Destructor 
}

//_______________________________________________________________________________
void AliMUONSimLoader::SetMUONData(AliMUONSimData * MUONData)
{
/// Set MUON data

  fMUONData = MUONData;
}
//_______________________________________________________________________________
AliMUONSimData * AliMUONSimLoader::GetMUONData()
{
/// Get MUON data

  return fMUONData;
}

