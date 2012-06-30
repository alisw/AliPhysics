// $Id$

/**************************************************************************
 * This file is property of and copyright by the                          * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTOfflineDataSource.cxx
/// @author Matthias Richter
/// @date   
/// @brief  AliRoot data source component base class.
///

#include "AliHLTOfflineDataSource.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOfflineDataSource)

AliHLTOfflineDataSource::AliHLTOfflineDataSource()
{
  // The class implements an AliRoot data source component base class.
  // Main purpose is to provide the access to RunLoader and RawReader
  // of the AliRoot simulation/reconstruction
  // Child classes implement the normal component functions to
  // facilitate the interface
  Register(this);
}

AliHLTOfflineDataSource::~AliHLTOfflineDataSource()
{
  // destructor
  Unregister(this);
}
