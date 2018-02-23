#ifndef AliAODFmdCluster_H
#define AliAODFmdCluster_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD FMD cluster class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TRef.h>

#include "AliAODCluster.h"
#include "AliAODTrack.h"

class AliAODFmdCluster : public AliAODCluster {

 public:
  
  AliAODFmdCluster();
  AliAODFmdCluster(Int_t id,
		   Int_t nLabel,
		   Int_t *label,
		   Double_t energy,
		   Double_t x[3],
		   Double_t pid[9],
		   Char_t ttype=kUndef,
		   AliAODVertex *prodVertex=NULL,
		   AliAODTrack *primTrack=NULL);

   AliAODFmdCluster(Int_t id,
		    Int_t nLabel,
		    Int_t *label,
		    Float_t energy,
		    Float_t x[3],
		    Float_t pid[9],
		    Char_t ttype=kUndef,
		    AliAODVertex *prodVertex=NULL,
		    AliAODTrack *primTrack=NULL);

  virtual ~AliAODFmdCluster();
  AliAODFmdCluster(const AliAODFmdCluster& trk); 
  AliAODFmdCluster& operator=(const AliAODFmdCluster& trk);

  AliAODVertex *GetProdVertex() const { return (AliAODVertex*)fProdVertex.GetObject(); }
  AliAODTrack *GetPrimTrack() const { return (AliAODTrack*)fPrimTrack.GetObject(); }
  
  void SetProdVertex(AliAODVertex *vertex) { fProdVertex = vertex; }
  void SetPrimTrack(AliAODTrack *ptrack) { fPrimTrack = ptrack; }


 private :

  TRef fProdVertex;     // vertex of origin
  TRef fPrimTrack;      // primary track associated with this cluster

  ClassDef(AliAODFmdCluster,1);
};

#endif
