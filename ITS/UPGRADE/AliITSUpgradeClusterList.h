#ifndef ITSUPGRADE_CLUSTER_LIST_H
#define ITSUPGRADE_CLUSTER_LIST_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//////////////////////////////////////////////////////////////////////
// Author: A.Mastroserio, C.Terrevoli                               //
//         annalisa.mastroserio@cern.ch			            //	
//         cristina.terrevoli@ba.infn.it                            //     
//							            //			
//   This is class implements the use of a list of clusters.        //
//                                                                  //
//////////////////////////////////////////////////////////////////////  

/* $Id$ */
#include <Rtypes.h>

class AliITSUpgradeClusterListNode;

class AliITSUpgradeClusterList {

 public:
  AliITSUpgradeClusterList();
 AliITSUpgradeClusterList(AliITSUpgradeClusterListNode* first,UInt_t nrEntries);
  AliITSUpgradeClusterList(const AliITSUpgradeClusterList& ilist);
  virtual ~AliITSUpgradeClusterList();
  AliITSUpgradeClusterList& operator=(const AliITSUpgradeClusterList& ilist);

  void                   Clear();
  AliITSUpgradeClusterList*     Clone() const;
 Bool_t                 Insert(Float_t col, Float_t row, UShort_t size, UShort_t widthZ, UShort_t widthPhi, UShort_t type, UShort_t charge, Int_t labels[12]);

  UInt_t   GetNrEntries() const {return fNrEntries;}
  Float_t  GetColIndex(UInt_t index);
  Float_t  GetRowIndex(UInt_t index);
  UShort_t GetSizeIndex(UInt_t index);
  UShort_t GetWidthZIndex(UInt_t index);
  UShort_t GetWidthPhiIndex(UInt_t index);
  UShort_t GetTypeIndex(UInt_t index);
  UShort_t GetCharge(UInt_t index);
  Int_t *  GetLabels(UInt_t index);

 private:
  UInt_t                 fNrEntries;
  AliITSUpgradeClusterListNode* fFirst;
  AliITSUpgradeClusterListNode* fLast;

  Bool_t                  fFastAccess;
  AliITSUpgradeClusterListNode** fFastAccessArray;
  UInt_t                  fDummyIndex;


  void                   ClearNode(AliITSUpgradeClusterListNode* &node); // delete this node and all after
  AliITSUpgradeClusterListNode* CloneNode(AliITSUpgradeClusterListNode* node) const;

  void                   ClearFastAccess();
  void                   InitFastAccess();
  void                   InitFastAccessNode(AliITSUpgradeClusterListNode* node);
  void                   SetLastNode();

};

#endif
