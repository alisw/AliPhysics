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

/*
$Log$
Revision 1.2.4.4  2000/06/26 07:39:42  kowal2
Changes to obey the coding rules

Revision 1.2.4.3  2000/06/25 08:38:41  kowal2
Splitted from AliTPCtracking

Revision 1.2.4.2  2000/06/14 16:45:13  kowal2
Improved algorithms. Compiler warnings removed.

Revision 1.2.4.1  2000/06/09 07:09:29  kowal2

Clustering and tracking classes are splitted from the simulation ones

Revision 1.2  2000/04/17 09:37:33  kowal2
removed obsolete AliTPCDigitsDisplay.C

Revision 1.1.4.2  2000/04/10 11:34:02  kowal2

Clusters handling in a new data structure

*/

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
#include "AliSegmentArray.h"

#include "AliCluster.h"
#include "AliClusters.h"
#include "AliDetectorParam.h"
#include "AliClustersArray.h"



ClassImp(AliClustersArray)
//

AliClustersArray::AliClustersArray()
{
  //
  //Default constructor
  //
  fParam = 0;
  fClusterType = 0;
}

Bool_t  AliClustersArray::SetClusterType(Text_t * classname) 
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

Bool_t AliClustersArray::Setup(const AliDetectorParam *param)
{
  //
  //make copy of param object
  
  return kTRUE;
}

Bool_t AliClustersArray::SetParam(AliDetectorParam * param)
{
  return kTRUE;
}

Bool_t AliClustersArray::SetFinder(AliClustersFinder * finder)
{
  return kTRUE;
}
