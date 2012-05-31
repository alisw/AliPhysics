#ifndef AliAODPmdCluster_H
#define AliAODPmdCluster_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD PMD cluster class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TRef.h>

#include "AliAODCluster.h"
class AliAODPmdCluster : public AliAODCluster {

 public:
  
  AliAODPmdCluster();
  AliAODPmdCluster(Int_t id,
		   Int_t nLabel,
		   Int_t *label,
		   Double_t energy,
		   Double_t x[3],
		   Double_t pid[13],
		   Char_t ttype=kUndef,
		   UInt_t selectInfo=0,
		   AliAODPmdCluster* assoc=NULL);
  
  AliAODPmdCluster(Int_t id,
		   Int_t nLabel,
		   Int_t *label,
		   Float_t energy,
		   Float_t x[3],
		   Float_t pid[13],
		   Char_t ttype=kUndef,
		   UInt_t selectInfo=0,
		   AliAODPmdCluster* assoc=NULL);
  
  virtual ~AliAODPmdCluster();
  AliAODPmdCluster(const AliAODPmdCluster& clus); 
  AliAODPmdCluster& operator=(const AliAODPmdCluster& clus);

  AliAODPmdCluster *GetAssocCluster() const { return (AliAODPmdCluster*)fAssocCluster.GetObject(); }
  void SetAssocCluster(AliAODPmdCluster *cluster) { fAssocCluster = cluster; }

 private :

  TRef fAssocCluster;   // cluster of other layer associated with this cluster

  ClassDef(AliAODPmdCluster,1);
};

#endif
