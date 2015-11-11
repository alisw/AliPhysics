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

/// \class AliTPCClustersRow
/// \brief Time Projection Chamber AliTPCClusterRow  objects
///
/// -   clusters for given segment of TPC
///
/// \author Marian Ivanov , GSI Darmstadt

#include <TClass.h>
#include "AliClusters.h"
#include "AliTPCclusterMI.h"
#include "AliTPCClustersRow.h"
#include "AliLog.h"
#include <TDirectory.h>
#include <TClonesArray.h>


const Int_t kDefSize = 1;  ///< defalut size


/// \cond CLASSIMP
ClassImp(AliTPCClustersRow)
/// \endcond


// code audit 2015-11-06 the class is in its current imlementation bound to
// AliTPCclusterMI as element class, corresponding code and checks have been
// added

//*****************************************************************************
//
//_____________________________________________________________________________
AliTPCClustersRow::AliTPCClustersRow() : AliClusters("AliTPCclusterMI")
{
  //
  //default constructor
  fNclusters=0;
}

//____________________________________________________________________________
AliTPCClustersRow::AliTPCClustersRow(const char *classname) : AliClusters("AliTPCclusterMI")
{
 /// special constructor
  TString cmpstr(classname);
  if (cmpstr.CompareTo("AliTPCclusterMI")) {
    AliFatal("Class AliTPCClustersRow is specifically bound to AliTPCclusterMI as element class");
  }

 fNclusters=0;
}

//_____________________________________________________________________________
TObject *AliTPCClustersRow::InsertCluster(const TObject *c)
{
  /// Add a simulated cluster copy to the list

  // code audit 2015-11-06 usage of name of fClass to create the TClonesArray
  // does not make sense, because the rest of the function is hardwired
  if (fClass==0) {
    Error("AliClusters", "class type not specified");
    return 0;
  }
  if(!fClusters) fClusters=new TClonesArray("AliTPCclusterMI",1000);
  TClonesArray &lclusters = *fClusters;
  return new(lclusters[fNclusters++]) AliTPCclusterMI(*((AliTPCclusterMI*)c));
}
//__________________________________________________________________________


TObject *AliTPCClustersRow::Append(){
 /// create new object return pointer to this object

 return fClusters->operator[](fClusters->GetEntriesFast());
}

