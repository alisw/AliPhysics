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

#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliAODHeader.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODCluster.h"
#include "AliAODJet.h"
#include "AliAODTracklets.h"

class TTree;

class AliAODEvent : public AliVEvent {

 public :
  enum AODListIndex_t {kAODHeader,
		       kAODTracks,
		       kAODVertices,
		       kAODClusters,
		       kAODJets,
		       kAODTracklets,
		       kAODListN
  };

  AliAODEvent();
  virtual ~AliAODEvent();

  //AliAODEvent(const AliAODEvent& aodevent);  // not implemented
  //AliAODEvent& operator=(const AliAODEvent& aodevent);  // not implemented

  void          AddObject(TObject *obj);
  TObject      *FindListObject(const char *objName) ;
  TList        *GetList()                const { return fAODObjects; }

  // -- Header
  AliAODHeader *GetHeader()              const { return fHeader; }
  void          AddHeader(const AliAODHeader* hdx)
    {
	delete fHeader; fHeader = new AliAODHeader(*hdx);
	(fAODObjects->FirstLink())->SetObject(fHeader);
    }

  // setters and getters for header information
  void     SetRunNumber(Int_t n) {if (fHeader) fHeader->SetRunNumber(n);}
  void     SetPeriodNumber(UInt_t n){if (fHeader) fHeader->SetPeriodNumber(n);}
  void     SetOrbitNumber(UInt_t n) {if (fHeader) fHeader->SetOrbitNumber(n);}
  void     SetBunchCrossNumber(UShort_t n) {if (fHeader) fHeader->SetBunchCrossNumber(n);}
  void     SetMagneticField(Double_t mf){if (fHeader) fHeader->SetMagneticField(mf);}

  Int_t    GetRunNumber() const {return fHeader ? fHeader->GetRunNumber() : -999;}
  UInt_t   GetPeriodNumber() const {return fHeader ? fHeader->GetPeriodNumber() : 0;}
  UInt_t   GetOrbitNumber() const {return fHeader ? fHeader->GetOrbitNumber() : 0;}
  UShort_t GetBunchCrossNumber() const {return fHeader ? fHeader->GetBunchCrossNumber() : 0;}
  Double_t GetMagneticField() const {return fHeader ? fHeader->GetMagneticField() : -999.;}
  
  void      SetEventType(UInt_t eventType){fHeader->SetEventType(eventType);}
  void      SetTriggerMask(ULong64_t n) {fHeader->SetTriggerMask(n);}
  void      SetTriggerCluster(UChar_t n) {fHeader->SetTriggerCluster(n);}

  UInt_t    GetEventType()  const { return fHeader ? fHeader->GetEventType() : 0;}
  ULong64_t GetTriggerMask() const {return fHeader ? fHeader->GetTriggerMask() : 0;}
  UChar_t   GetTriggerCluster() const {return fHeader ? fHeader->GetTriggerCluster() : 0;}

  Double_t  GetZDCN1Energy()        const { return fHeader ? fHeader->GetZDCN1Energy() : -999.; }
  Double_t  GetZDCP1Energy()        const { return fHeader ? fHeader->GetZDCP1Energy() : -999.; }
  Double_t  GetZDCN2Energy()        const { return fHeader ? fHeader->GetZDCN2Energy() : -999.; }
  Double_t  GetZDCP2Energy()        const { return fHeader ? fHeader->GetZDCP2Energy() : -999.; }
  Double_t  GetZDCEMEnergy()        const { return fHeader ? fHeader->GetZDCEMEnergy() : -999.; }

  // -- Tracks
  TClonesArray *GetTracks()              const { return fTracks; }
  Int_t         GetNTracks()             const { return fTracks->GetEntriesFast(); }
  Int_t         GetNumberOfTracks()      const { return GetNTracks(); }
  AliAODTrack  *GetTrack(Int_t nTrack)   const { return (AliAODTrack*)fTracks->UncheckedAt(nTrack); }
  Int_t         AddTrack(const AliAODTrack* trk)
    {new((*fTracks)[fTracks->GetEntriesFast()]) AliAODTrack(*trk); return fTracks->GetEntriesFast()-1;}
  Int_t         GetMuonTracks(TRefArray *muonTracks) const;

  // -- Vertex
  TClonesArray *GetVertices()            const { return fVertices; }
  Int_t         GetNVertices()           const { return fVertices->GetEntriesFast(); }
  AliAODVertex *GetVertex(Int_t nVertex) const { return (AliAODVertex*)fVertices->UncheckedAt(nVertex); }
  Int_t         AddVertex(const AliAODVertex* vtx)
  {new((*fVertices)[fVertices->GetEntriesFast()]) AliAODVertex(*vtx); return fVertices->GetEntriesFast()-1;}
  virtual AliAODVertex *GetPrimaryVertex() const { return GetVertex(0); }
  

  // -- Cluster
  TClonesArray *GetClusters()            const { return fClusters; }
  Int_t         GetNClusters()           const { return fClusters->GetEntriesFast(); }
  AliAODCluster *GetCluster(Int_t nCluster) const { return (AliAODCluster*)fClusters->UncheckedAt(nCluster); }
  Int_t         AddCluster(const AliAODCluster* vtx)
  {new((*fClusters)[fClusters->GetEntriesFast()]) AliAODCluster(*vtx); return fClusters->GetEntriesFast()-1;}

  // -- Jet
  TClonesArray *GetJets()            const { return fJets; }
  Int_t         GetNJets()           const { return fJets->GetEntriesFast(); }
  AliAODJet    *GetJet(Int_t nJet) const { return (AliAODJet*)fJets->UncheckedAt(nJet); }
  Int_t         AddJet(const AliAODJet* vtx)
    {new((*fJets)[fJets->GetEntriesFast()]) AliAODJet(*vtx); return fJets->GetEntriesFast()-1;}

  // -- Tracklets
  AliAODTracklets *GetTracklets() const { return fTracklets; }

  // -- Services
  void    CreateStdContent();
  void    SetStdNames();
  void    GetStdContent();
  void    ResetStd(Int_t trkArrSize = 0, Int_t vtxArrSize = 0);
  void    ClearStd();
  void    ReadFromTree(TTree *tree);
  const void WriteToTree(TTree* tree) const {tree->Branch(fAODObjects);}

  void  Print(Option_t *option="") const;

 private :

  TList *fAODObjects; //  list of AODObjects
  
  // standard content
  AliAODHeader    *fHeader;    //! event information
  TClonesArray    *fTracks;    //! charged tracks
  TClonesArray    *fVertices;  //! vertices
  TClonesArray    *fClusters;  //! neutral particles
  TClonesArray    *fJets;      //! jets
  AliAODTracklets *fTracklets; //! SPD tracklets

  static const char* fAODListName[kAODListN]; //!

  ClassDef(AliAODEvent,3);
};

#endif
