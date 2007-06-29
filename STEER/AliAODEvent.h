#ifndef AliAODEvent_H
#define AliAODEvent_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TBuffer.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TTree.h>
#include <TNamed.h>

#include "AliAODHeader.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODCluster.h"
#include "AliAODJet.h"

class TTree;

class AliAODEvent : public TObject {

 public :
  AliAODEvent();
  virtual ~AliAODEvent();

  //AliAODEvent(const AliAODEvent& aodevent);  // not implemented
  //AliAODEvent& operator=(const AliAODEvent& aodevent);  // not implemented

  void          AddObject(TObject *obj);
  TObject      *GetObject(const char *objName) const;
  TList        *GetList()                const { return fAODObjects; }

  // -- Header
  AliAODHeader *GetHeader()              const { return (AliAODHeader*)fHeader; }
  void          AddHeader(const AliAODHeader* hdx)
    {
	delete fHeader; fHeader = new AliAODHeader(*hdx);
	(fAODObjects->FirstLink())->SetObject(fHeader);
    }

  // -- Tracks
  TClonesArray *GetTracks()              const { return fTracks; }
  Int_t         GetNTracks()             const { return fTracks->GetEntriesFast(); }
  AliAODTrack  *GetTrack(Int_t nTrack)   const { return (AliAODTrack*)fTracks->At(nTrack); }
  void          AddTrack(const AliAODTrack* trk)
    {new((*fTracks)[fTracks->GetEntries()]) AliAODTrack(*trk);}
  Int_t         GetMuonTracks(TRefArray *muonTracks) const;

  // -- Vertex
  TClonesArray *GetVertices()            const { return fVertices; }
  Int_t         GetNVertices()           const { return fVertices->GetEntriesFast(); }
  AliAODVertex *GetVertex(Int_t nVertex) const { return (AliAODVertex*)fVertices->At(nVertex); }
  void          AddVertex(const AliAODVertex* vtx)
    {new((*fVertices)[fVertices->GetEntries()]) AliAODVertex(*vtx);}
  virtual AliAODVertex *GetPrimaryVertex() const { return GetVertex(0); }
  

  // -- Cluster
  TClonesArray *GetClusters()            const { return fClusters; }
  Int_t         GetNClusters()           const { return fClusters->GetEntriesFast(); }
  AliAODCluster *GetCluster(Int_t nCluster) const { return (AliAODCluster*)fClusters->At(nCluster); }
  void          AddCluster(const AliAODCluster* vtx)
    {new((*fClusters)[fClusters->GetEntries()]) AliAODCluster(*vtx);}

  // -- Jet
  TClonesArray *GetJets()            const { return fJets; }
  Int_t         GetNJets()           const { return fJets->GetEntriesFast(); }
  AliAODJet *GetJet(Int_t nJet) const { return (AliAODJet*)fJets->At(nJet); }
  void          AddJet(const AliAODJet* vtx)
    {new((*fJets)[fJets->GetEntries()]) AliAODJet(*vtx);}

  // -- Services
  void    CreateStdContent();
  void    GetStdContent();
  void    ResetStd(Int_t trkArrSize = 0, Int_t vtxArrSize = 0);
  void    ClearStd();
  void    ReadFromTree(TTree *tree);
  void    WriteToTree(TTree* tree) {tree->Branch(fAODObjects);}
 private :

  AliAODEvent(const AliAODEvent&); // Not implemented
  AliAODEvent& operator=(const AliAODEvent&); // Not implemented

  TList *fAODObjects; //  list of AODObjects
  
  // standard content
  AliAODHeader  *fHeader;   //! event information
  TClonesArray  *fTracks;   //! charged tracks
  TClonesArray  *fVertices; //! vertices
  TClonesArray  *fClusters; //! neutral particles
  TClonesArray  *fJets;     //! jets

  ClassDef(AliAODEvent,1);
};

#endif
