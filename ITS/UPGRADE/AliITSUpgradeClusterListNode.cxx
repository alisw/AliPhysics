/**************************************************************************
 * Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////
// Author:A.Mastroserio C.Terrevoli                                 //
//        annalisa.mastroserio@cern.ch                              //
//        cristina.terrevoli@ba.infn.it                             //
//                                                                  //
// Implementation of the nodes to put in the cluster list           //
//                              			            //
//////////////////////////////////////////////////////////////////////   

#include "AliITSUpgradeClusterListNode.h"

AliITSUpgradeClusterListNode::AliITSUpgradeClusterListNode():
  fCol(0),
  fRow(0),
  fSize(0),
  fWidthZ(0),
  fWidthPhi(0),
  fType(0),
  fCharge(0),
  fLastDigitLabel(0),
  fNext(NULL)
{
for(Int_t i=0; i < kMaxLab; i++) fDigitLabel[i]=-2;
}

AliITSUpgradeClusterListNode::AliITSUpgradeClusterListNode(Float_t col, Float_t row, UShort_t size, UShort_t widthZ, UShort_t widthPhi, UShort_t type, UShort_t charge, AliITSUpgradeClusterListNode* next):
  fCol(col),
  fRow(row),
  fSize(size),
  fWidthZ(widthZ),
  fWidthPhi(widthPhi),
  fType(type),
  fCharge(charge),
  fLastDigitLabel(0),
  fNext(next)
{
for(Int_t i=0; i < kMaxLab; i++) fDigitLabel[i]=-2;

}

AliITSUpgradeClusterListNode::AliITSUpgradeClusterListNode(const AliITSUpgradeClusterListNode& obj):
  fCol(obj.fCol),
  fRow(obj.fRow),
  fSize(obj.fSize),
  fWidthZ(obj.fWidthZ),
  fWidthPhi(obj.fWidthPhi),
  fType(obj.fType),
  fCharge(obj.fCharge),
  fLastDigitLabel(obj.fLastDigitLabel),
  fNext(obj.fNext)
{
  // copy constructor
for(Int_t i=0; i< kMaxLab; i++) fDigitLabel[i]=obj.fDigitLabel[i];
}

AliITSUpgradeClusterListNode::~AliITSUpgradeClusterListNode() 
{}

AliITSUpgradeClusterListNode& AliITSUpgradeClusterListNode::operator=(const AliITSUpgradeClusterListNode& obj) 
{
  // assignment operator
  if (this!=&obj) {
    fCol = obj.fCol;
    fRow = obj.fRow;
    fSize = obj.fSize;
    fWidthZ = obj.fWidthZ;
    fWidthPhi = obj.fWidthPhi;
    fType = obj.fType;
    fCharge = obj.fCharge;
    fNext = obj.fNext;
  }
  return *this;
}

