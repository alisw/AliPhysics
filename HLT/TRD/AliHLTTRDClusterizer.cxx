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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// HLT TRD cluster finder                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliHLTTRDClusterizer.h"
#include "AliTRDgeometry.h"
#include "AliTRDcluster.h"
#include "AliTRDReconstructor.h"
#include <TClonesArray.h>

ClassImp(AliHLTTRDClusterizer);

//_____________________________________________________________________________
AliHLTTRDClusterizer::AliHLTTRDClusterizer(const AliTRDReconstructor *const rec)
  :AliTRDclusterizer(rec)
  , fMemBlock(NULL)
{
  //
  // AliHLTTRDClusterizer default constructor
  //
}

//_____________________________________________________________________________
AliHLTTRDClusterizer::AliHLTTRDClusterizer(const Text_t *const name, const Text_t *const title, const AliTRDReconstructor *const rec)
  : AliTRDclusterizer(name,title,rec)
  , fMemBlock(NULL)
{
  //
  // AliHLTTRDClusterizer constructor
  //
}

//_____________________________________________________________________________
AliHLTTRDClusterizer::AliHLTTRDClusterizer(const AliHLTTRDClusterizer& c)
  : AliTRDclusterizer(c)
  , fMemBlock(NULL)
{
  //
  // AliHLTTRDClusterizer copy constructor
  //
}

//_____________________________________________________________________________
AliHLTTRDClusterizer& AliHLTTRDClusterizer::operator=(const AliHLTTRDClusterizer& c)
{
  //
  // Assignment operator
  //
  
  if(this!=&c) 
    c.Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void AliHLTTRDClusterizer::Copy(TObject& c) const
{
  //
  // Copy function
  //

  ((AliHLTTRDClusterizer&)c).fMemBlock  = NULL;
}

//_____________________________________________________________________________
void AliHLTTRDClusterizer::AddClusterToArray(AliTRDcluster *cluster)
{
  //
  // Add a cluster to the array
  //

  AliHLTTRDCluster *ptr = &(((AliHLTTRDCluster*)GetMemBlock())[fNoOfClusters]);
  new(ptr) AliHLTTRDCluster(cluster);
}
