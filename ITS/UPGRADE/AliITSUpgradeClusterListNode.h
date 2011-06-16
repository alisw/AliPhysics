#ifndef ITSUPGRADE_CLUSTER_LISTNODE_H
#define ITSUPGRADE_CLUSTER_LISTNODE_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
////////////////////////////////////////////////////////////////////////
// Author:A.Mastroserio C.Terrevoli 				      //	
//        annalisa.mastroserio@cern.ch  			      //	
//        cristina.terrevoli@ba.infn.it		                      //
// 								      //	
//  Class for the nodes to put in the cluster list		      //
//                                                                    //
////////////////////////////////////////////////////////////////////////  

/* $Id$ */

#include <Rtypes.h>

class AliITSUpgradeClusterListNode {

 public:
 AliITSUpgradeClusterListNode();
  AliITSUpgradeClusterListNode(Float_t col, Float_t row, UShort_t size, UShort_t widthZ, UShort_t widthPhi, UShort_t type,UShort_t charge,AliITSUpgradeClusterListNode* next);
  AliITSUpgradeClusterListNode(const AliITSUpgradeClusterListNode& obj);
  virtual ~AliITSUpgradeClusterListNode();
  AliITSUpgradeClusterListNode& operator=(const AliITSUpgradeClusterListNode& obj);

  Float_t  Col() const {return fCol;}
  Float_t  Row() const {return fRow;}
  UShort_t Size() const {return fSize;}
  UShort_t WidthZ() const {return fWidthZ;}
  UShort_t WidthPhi() const {return fWidthPhi;}
  UShort_t Type() const {return fType;}
  UShort_t Charge() const {return fCharge;}
  AliITSUpgradeClusterListNode*& Next() {return fNext;}

  void       SetCol(Float_t val) {fCol=val;}
  void       SetRow(Float_t val) {fRow=val;}
  void       SetSize(UShort_t val) {fSize=val;}
  void       SetWidthZ(UShort_t val) {fWidthZ=val;}
  void       SetWidthPhi(UShort_t val) {fWidthPhi=val;}
  void       SetType(UShort_t val) {fType=val;}
  void       SetNext(AliITSUpgradeClusterListNode* obj) {fNext = obj;}
 
  void       AddDigitLabel(Int_t label){fDigitLabel[fLastDigitLabel]=label; fLastDigitLabel++;}
  Int_t      GetLastIndex() {return fLastDigitLabel;}
  Int_t*     GetLabels() {return fDigitLabel;}



 private:
  Float_t  fCol;
  Float_t  fRow;
  UShort_t fSize;
  UShort_t fWidthZ;
  UShort_t fWidthPhi;
  UShort_t fType;
  UShort_t fCharge;
  Int_t  fLastDigitLabel;   // last meaningful position in fDigitLabel
  AliITSUpgradeClusterListNode* fNext;

  enum {kMaxLab=24*12}; // maximum number of MC labels associated to the cluster
  Int_t  fDigitLabel[kMaxLab]; // needed to attach MC truth to the cluster
};

#endif

