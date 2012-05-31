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
//  AliClustersArray  object                                                 //
//
//  Origin: Marian Ivanov , GSI Darmstadt                                    //
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TClass.h"
#include  <TROOT.h>
#include "AliSegmentID.h"
#include "TObjArray.h"

#include "AliClusters.h"
#include "AliDetectorParam.h"
#include "AliClustersArray.h"



ClassImp(AliClustersArray)
//

  AliClustersArray::AliClustersArray():AliSegmentArray(),
				       fParam(0),
				       fClFinder(0),
				       fClusterType(0)
{
  //
  //Default constructor
  //

}

Bool_t  AliClustersArray::SetClusterType(const char * classname) 
{
  //
  //set type of Clusters
  //
  if ( fClusterType !=0 ) {
    delete fClusterType;
    fClusterType = 0;
  }

  if (!gROOT)
    ::Fatal("AliClustersArray", "ROOT system not initialized");
   
   fClusterType = gROOT->GetClass(classname);
   if (!fClusterType) {
      Error("AliClustersArray", "%s is not a valid class name", classname);
      return kFALSE;
   }
   if (!fClusterType->InheritsFrom(TObject::Class())) {
      Error("AliClustersArray", "%s does not inherit from TObject", classname);
      return kFALSE;
   }     
  return kTRUE;
}

Bool_t AliClustersArray::Setup(const AliDetectorParam */*param*/)
{
  //
  //make copy of param object
  
  return kTRUE;
}

Bool_t AliClustersArray::SetParam(AliDetectorParam * /*param*/)
{
  return kTRUE;
}

Bool_t AliClustersArray::SetFinder(AliClustersFinder * /*finder*/)
{
  return kTRUE;
}
