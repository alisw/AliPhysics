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
#include "AliLog.h"

ClassImp(AliTPCclusterMI)


AliTPCclusterMI::AliTPCclusterMI(Bool_t withInfo):
  AliCluster(),
  fX(0),
  fQ(0),
  fType(0),
  fMax(0),
  fUsed(0),
  fDetector(0),
  fRow(0),
  fTimeBin(0),
  fPad(0),
  fInfo(0)
{
  //
  // default constructor
  //
  if (withInfo) fInfo = new AliTPCclusterInfo;
}

AliTPCclusterMI::AliTPCclusterMI(const AliTPCclusterMI & cluster):
  AliCluster(cluster),
  fX(cluster.fX),
  fQ(cluster.fQ),
  fType(cluster.fType),
  fMax(cluster.fMax),
  fUsed(cluster.fUsed),
  fDetector(cluster.fDetector),
  fRow(cluster.fRow),
  fTimeBin(cluster.fTimeBin),
  fPad(cluster.fPad),
  fInfo(0)
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

  (AliCluster&)(*this) = (AliCluster&)cluster;
  fX    = cluster.fX;
  fQ    = cluster.fQ;
  fType = cluster.fType;
  fMax  = cluster.fMax;
  fUsed = cluster.fUsed;
  fDetector = cluster.fDetector;
  fRow  = cluster.fRow;
  fTimeBin = cluster.fTimeBin;
  fPad     = cluster.fPad;
  fInfo = 0;
  if (cluster.fInfo) fInfo = new AliTPCclusterInfo(*(cluster.fInfo));
}




AliTPCclusterMI::AliTPCclusterMI(Int_t *lab, Float_t *hit) : 
  AliCluster(lab,hit),
  fX(0),
  fQ(0),
  fType(0),
  fMax(0),
  fUsed(0),
  fDetector(0),
  fRow(0),
  fInfo(0)
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
  return (o2->GetY()>fY)? -1:1; 
}
