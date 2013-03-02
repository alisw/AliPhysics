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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber AliTPCClusterRow  objects
//  -   clusters for given segment of TPC                                //
//
//  Origin: Marian Ivanov , GSI Darmstadt
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
#include <TClass.h>
#include "AliClusters.h"
#include "AliTPCclusterMI.h"
#include "AliTPCClustersRow.h"
#include <TDirectory.h>
#include <TClonesArray.h>


const Int_t kDefSize = 1;  //defalut size


ClassImp(AliTPCClustersRow) 


//*****************************************************************************
//
//_____________________________________________________________________________
AliTPCClustersRow::AliTPCClustersRow() 
{  
  //
  //default constructor
  fNclusters=0;
}

//____________________________________________________________________________
AliTPCClustersRow::AliTPCClustersRow(const char *classname) : AliClusters(classname)
{
 // special constructor
 fNclusters=0;
}

//_____________________________________________________________________________
TObject *AliTPCClustersRow::InsertCluster(const TObject *c) 
{    
  //
  // Add a simulated cluster copy to the list
  //
  if (fClass==0) {
    Error("AliClusters", "class type not specified");
    return 0;
  }
  if(!fClusters) fClusters=new TClonesArray(fClass->GetName(),1000);
  TClonesArray &lclusters = *fClusters;
  return new(lclusters[fNclusters++]) AliTPCclusterMI(*((AliTPCclusterMI*)c));
}
//__________________________________________________________________________


TObject *AliTPCClustersRow::Append(){
 //create new object return pointer to this object
 return fClusters->operator[](fClusters->GetEntriesFast());
}

