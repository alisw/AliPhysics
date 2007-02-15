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
#include <TNamed.h>

#include "AliAODHeader.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODNeutral.h"
#include "AliAODJet.h"

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
    {delete fHeader; fHeader=new AliAODHeader(*hdx);}

  // -- Tracks
  TClonesArray *GetTracks()              const { return fTracks; }
  Int_t         GetNTracks()             const { return fTracks->GetEntriesFast(); }
  AliAODTrack  *GetTrack(Int_t nTrack)   const { return (AliAODTrack*)fTracks->At(nTrack); }
  void          AddTrack(const AliAODTrack* trk)
    {new((*fTracks)[fTracks->GetEntries()]) AliAODTrack(*trk);}

  // -- Vertex
  TClonesArray *GetVertices()            const { return fVertices; }
  Int_t         GetNVertices()           const { return fVertices->GetEntriesFast(); }
  AliAODVertex *GetVertex(Int_t nVertex) const { return (AliAODVertex*)fVertices->At(nVertex); }
  void          AddVertex(const AliAODVertex* vtx)
    {new((*fVertices)[fVertices->GetEntries()]) AliAODVertex(*vtx);}

  // -- Neutral
  TClonesArray *GetNeutrals()            const { return fNeutrals; }
  Int_t         GetNNeutrals()           const { return fNeutrals->GetEntriesFast(); }
  AliAODNeutral *GetNeutral(Int_t nNeutral) const { return (AliAODNeutral*)fNeutrals->At(nNeutral); }
  void          AddNeutral(const AliAODNeutral* vtx)
    {new((*fNeutrals)[fNeutrals->GetEntries()]) AliAODNeutral(*vtx);}

  // -- Jet
  TClonesArray *GetJets()            const { return fJets; }
  Int_t         GetNJets()           const { return fJets->GetEntriesFast(); }
  AliAODJet *GetJet(Int_t nJet) const { return (AliAODJet*)fJets->At(nJet); }
  void          AddJet(const AliAODJet* vtx)
    {new((*fJets)[fJets->GetEntries()]) AliAODJet(*vtx);}

  // -- Services
  void CreateStdContent();
  void GetStdContent() const;
  void ResetStd(Int_t trkArrSize = 0, Int_t vtxArrSize = 0);

 private :

  AliAODEvent(const AliAODEvent&); // Not implemented
  AliAODEvent& operator=(const AliAODEvent&); // Not implemented

  TList *fAODObjects; // list of AODObjects

  // standard content
  mutable AliAODHeader  *fHeader;   //! event information
  mutable TClonesArray  *fTracks;   //! charged tracks
  mutable TClonesArray  *fVertices; //! vertices
  mutable TClonesArray  *fNeutrals; //! neutral particles
  mutable TClonesArray  *fJets;     //! jets

  ClassDef(AliAODEvent,1);
};

#endif
