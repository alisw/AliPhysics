// -*- mode: C++ -*- 
#ifndef ALIVEVENT_H
#define ALIVEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliVEvent
//      
// Origin: Markus Oldenburg, CERN, Markus.Oldenburg@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TTree.h>
#include "AliVHeader.h"
#include "AliVParticle.h"
#include "AliVVertex.h"

class AliVEvent : public TObject {

public:

  AliVEvent() { }
  virtual ~AliVEvent() { } 
  AliVEvent(const AliVEvent& vEvnt); 
  AliVEvent& operator=(const AliVEvent& vEvnt);

  // Services
  virtual void AddObject(TObject* obj) = 0;
  virtual TObject* FindListObject(const char *name) const = 0;
  virtual TList* GetList() const = 0;

  virtual void CreateStdContent() = 0;
  virtual void GetStdContent() = 0;

  virtual void ReadFromTree(TTree *tree, Option_t* opt) = 0;
  virtual void WriteToTree(TTree* tree) const = 0;

  virtual void Reset() = 0;
  //virtual void ResetStdContent() = 0;
  virtual void SetStdNames() = 0;

  virtual void Print(Option_t *option="") const = 0;

  // Header
  virtual AliVHeader* GetHeader() const = 0;

  // Delegated methods for fESDRun or AODHeader
  
  virtual void     SetRunNumber(Int_t n) = 0;
  virtual void     SetPeriodNumber(UInt_t n) = 0;
  virtual void     SetMagneticField(Double_t mf) = 0;
  
  virtual Int_t    GetRunNumber() const = 0;
  virtual UInt_t   GetPeriodNumber() const = 0;
  virtual Double_t GetMagneticField() const = 0;

  virtual Double_t GetDiamondX() const {return -999.;}
  virtual Double_t GetDiamondY() const {return -999.;}
  virtual void     GetDiamondCovXY(Float_t cov[3]) const
             {cov[0]=-999.; return;}

  // Delegated methods for fHeader
  virtual void      SetOrbitNumber(UInt_t n) = 0;
  virtual void      SetBunchCrossNumber(UShort_t n) = 0;
  virtual void      SetEventType(UInt_t eventType)= 0;
  virtual void      SetTriggerMask(ULong64_t n) = 0;
  virtual void      SetTriggerCluster(UChar_t n) = 0;

  virtual UInt_t    GetOrbitNumber() const = 0;
  virtual UShort_t  GetBunchCrossNumber() const = 0;
  virtual UInt_t    GetEventType()  const = 0;
  virtual ULong64_t GetTriggerMask() const = 0;
  virtual UChar_t   GetTriggerCluster() const = 0;

  virtual Double_t  GetZDCN1Energy() const = 0;
  virtual Double_t  GetZDCP1Energy() const = 0;
  virtual Double_t  GetZDCN2Energy() const = 0;
  virtual Double_t  GetZDCP2Energy() const = 0;
  virtual Double_t  GetZDCEMEnergy(Int_t i) const = 0;
 
  // Tracks
  virtual AliVParticle *GetTrack(Int_t i) const = 0;
  //virtual Int_t        AddTrack(const AliVParticle *t) = 0;
  virtual Int_t        GetNumberOfTracks() const = 0;
  virtual Int_t        GetNumberOfV0s() const = 0;
  virtual Int_t        GetNumberOfCascades() const = 0;

  // Primary vertex
  virtual const AliVVertex   *GetPrimaryVertex() const {return 0x0;}
  virtual Bool_t IsPileupFromSPD(Int_t /*minContributors*/, 
				 Double_t /*minZdist*/, 
				 Double_t /*nSigmaZdist*/, 
				 Double_t /*nSigmaDiamXY*/, 
				 Double_t /*nSigmaDiamZ*/)
				 const{
    return kFALSE;
  }
  //---------- end of new stuff


  /*  to be considered to go in here be implemented

  void SetPrimaryVertex(const AliESDVertex *vertex) {
    *fPrimaryVertex = *vertex;
    fPrimaryVertex->SetName("PrimaryVertex");// error prone use class wide names?
  }

  void SetMultiplicity(const AliMultiplicity *mul) {
    *fSPDMult = *mul;
    // CKB 
    //     new (&fSPDMult) AliMultiplicity(*mul);
  }
  const AliMultiplicity *GetMultiplicity() const {return fSPDMult;}
  
  
  AliESDMuonTrack *GetMuonTrack(Int_t i) const {
    return (AliESDMuonTrack *)fMuonTracks->UncheckedAt(i);
  }
  void AddMuonTrack(const AliESDMuonTrack *t) {
    TClonesArray &fmu = *fMuonTracks;
    new(fmu[fMuonTracks->GetEntriesFast()]) AliESDMuonTrack(*t);
  }

  AliESDv0 *GetV0(Int_t i) const {
    return (AliESDv0*)fV0s->UncheckedAt(i);
  }
  Int_t AddV0(const AliESDv0 *v);

  AliESDcascade *GetCascade(Int_t i) const {
    return (AliESDcascade *)fCascades->UncheckedAt(i);
  }
  void AddCascade(const AliESDcascade *c) {
    TClonesArray &fc = *fCascades;
    new(fc[fCascades->GetEntriesFast()]) AliESDcascade(*c);
  }

  AliESDkink *GetKink(Int_t i) const {
    return (AliESDkink *)fKinks->UncheckedAt(i);
  }
  Int_t AddKink(const AliESDkink *c);

  AliESDCaloCluster *GetCaloCluster(Int_t i) const {
    return (AliESDCaloCluster *)fCaloClusters->UncheckedAt(i);
  }
  Int_t AddCaloCluster(const AliESDCaloCluster *c);

  Int_t GetNumberOfMuonTracks() const {return fMuonTracks->GetEntriesFast();}
  Int_t GetNumberOfV0s()      const {return fV0s->GetEntriesFast();}
  Int_t GetNumberOfCascades() const {return fCascades->GetEntriesFast();}
  Int_t GetNumberOfKinks() const {return fKinks->GetEntriesFast();}
  Int_t GetNumberOfCaloClusters() const {return fCaloClusters->GetEntriesFast();}

  */

  ClassDef(AliVEvent,1)  // base class for AliEvent data
};
#endif 

