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
#include <TGeoMatrix.h>
#include "AliVHeader.h"
#include "AliVParticle.h"
#include "AliVVertex.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "TRefArray.h"
class AliCentrality;
class AliEventplane;
class AliVVZERO;

class AliVEvent : public TObject {

public:
  enum EOfflineTriggerTypes { 
      kMB           = BIT(0), // Minimum bias trigger, i.e. interaction trigger, offline SPD or V0 selection
                              // empty slot
      kMUON         = BIT(2), // Muon trigger, offline SPD or V0 selection
      kHighMult     = BIT(3), // High-multiplicity trigger (threshold defined online), offline SPD or V0 selection
      kEMC1         = BIT(4), // EMCAL trigger
      kUserDefined  = BIT(31), // Set when custom trigger classes are set in AliPhysicsSelection, offline SPD or V0 selection
      kAny          = 0xffffffff // to accept any trigger
  };

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

  // Calorimeter Clusters/Cells
  virtual AliVCluster *GetCaloCluster(Int_t)   const {return 0;}
  virtual Int_t GetNumberOfCaloClusters()      const {return 0;}
  virtual Int_t GetEMCALClusters(TRefArray *)  const {return 0;}
  virtual Int_t GetPHOSClusters (TRefArray *)  const {return 0;}
  virtual AliVCaloCells *GetEMCALCells()       const {return 0;}
  virtual AliVCaloCells *GetPHOSCells()        const {return 0;}
  const TGeoHMatrix* GetPHOSMatrix(Int_t /*i*/)    const {return NULL;}
  const TGeoHMatrix* GetEMCALMatrix(Int_t /*i*/)   const {return NULL;}

	
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

  virtual Bool_t IsPileupFromSPDInMultBins() const {
    return kFALSE;    
  }
  virtual AliCentrality* GetCentrality()                          = 0;
  virtual AliEventplane* GetEventplane()                          = 0;
  virtual Int_t        EventIndex(Int_t itrack)             const = 0;
  virtual Int_t        EventIndexForCaloCluster(Int_t iclu) const = 0;
  virtual Int_t        EventIndexForPHOSCell(Int_t icell)   const = 0;
  virtual Int_t        EventIndexForEMCALCell(Int_t icell)  const = 0;  

  virtual AliVVZERO *GetVZEROData() const = 0;

  ClassDef(AliVEvent,1)  // base class for AliEvent data
};
#endif 

