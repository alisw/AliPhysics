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

//-------------------------------------------------------
//          Implementation of the TPC cluser
//
//   Origin: Marian Ivanov   Marian.Ivanov@cern.ch
// 
//  AliTPC parallel tracker - 
//  Description of this class together with its intended usage
//  will follow shortly
//  
//-------------------------------------------------------

/* $Id$ */

#include "AliTPCclusterMI.h"
#include "AliTPCclusterInfo.h"
#include "AliGeomManager.h"
#include "AliLog.h"

ClassImp(AliTPCclusterMI)


AliTPCclusterMI::AliTPCclusterMI():
  AliCluster(),
  fInfo(0),
  fTimeBin(0),  //time bin coordinate
  fPad(0),  //pad coordinate
  fQ(0),       //Q of cluster (in ADC counts)  
  fMax(0),      //maximal amplitude in cluster
  fType(0),     //type of the cluster 0 means golden 
  fUsed(0),     //counter of usage  
  fDetector(0), //detector  number
  fRow(0)      //row number number
{
  //
  // default constructor
  //
}

AliTPCclusterMI::AliTPCclusterMI(const AliTPCclusterMI & cluster):
  AliCluster(cluster),
  fInfo(0),
  fTimeBin(cluster.fTimeBin),
  fPad(cluster.fPad),
  fQ(cluster.fQ),
  fMax(cluster.fMax),
  fType(cluster.fType),
  fUsed(cluster.fUsed),
  fDetector(cluster.fDetector),
  fRow(cluster.fRow)
{
  //
  // copy constructor
  // 
  //  AliInfo("Copy constructor\n");
  if (cluster.fInfo) fInfo = new AliTPCclusterInfo(*(cluster.fInfo));
}

AliTPCclusterMI & AliTPCclusterMI::operator = (const AliTPCclusterMI & cluster)
{
  //
  // assignment operator
  // 
  //  AliInfo("Asignment operator\n");
  if (this == &cluster) return (*this);

  (AliCluster&)(*this) = (AliCluster&)cluster;
  fQ    = cluster.fQ;
  fType = cluster.fType;
  fMax  = cluster.fMax;
  fUsed = cluster.fUsed;
  fDetector = cluster.fDetector;
  fRow  = cluster.fRow;
  fTimeBin = cluster.fTimeBin;
  fPad     = cluster.fPad;
  delete fInfo;
  fInfo = 0;
  if (cluster.fInfo) fInfo = new AliTPCclusterInfo(*(cluster.fInfo));
  return *this;
}




AliTPCclusterMI::AliTPCclusterMI(Int_t *lab, Float_t *hit) : 
  AliCluster(0,hit,0.,0.,lab),
  fInfo(0),
  fTimeBin(0),  //time bin coordinate
  fPad(0),  //pad coordinate
  fQ(0),       //Q of cluster (in ADC counts)  
  fMax(0),      //maximal amplitude in cluster
  fType(0),     //type of the cluster 0 means golden 
  fUsed(0),     //counter of usage  
  fDetector(0), //detector  number
  fRow(0)      //row number number
{
  //
  // constructor
  //
  fQ = (UShort_t)hit[4];
  fInfo = 0;
}

AliTPCclusterMI::~AliTPCclusterMI() {
  //
  // destructor
  //
  if (fInfo) delete fInfo;
  fInfo = 0;
}



Bool_t AliTPCclusterMI::IsSortable() const
{
  //
  //
  return kTRUE;

}

Int_t AliTPCclusterMI::Compare(const TObject* obj) const
{
  //
  // compare according y
  AliTPCclusterMI * o2 = (AliTPCclusterMI*)obj;
  return (o2->GetY()>GetY())? -1:1; 
}


void AliTPCclusterMI::SetDetector(Int_t detector){
  //
  // set volume ID 
  //  
  fDetector = (UChar_t)(detector%72);
  AliGeomManager::ELayerID id = (fDetector<36) ? 
    AliGeomManager::kTPC1 :AliGeomManager::kTPC2 ;
  Int_t modId = (fDetector<36)?fDetector: fDetector-36;
  SetVolumeId(AliGeomManager::LayerToVolUID(id,modId));  
}


void AliTPCclusterMI::SetInfo(AliTPCclusterInfo * info) {
  //
  //
  //
  if (fInfo) delete fInfo;
  fInfo = info;
}
