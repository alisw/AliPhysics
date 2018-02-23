/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD class to store tracklets
//     Author: Jan Fiete Grosse-Oetringhaus, CERN
//     Class created from AliMultiplicity
//-------------------------------------------------------------------------

#ifndef ALIAODTRACKLETS_H
#define ALIAODTRACKLETS_H

#include "AliVMultiplicity.h"
#include <TBits.h>

class AliAODTracklets : public AliVMultiplicity
{
 public:
  AliAODTracklets();
  AliAODTracklets(const char* name, const char* title);
  AliAODTracklets(const AliAODTracklets& evt); 
  AliAODTracklets& operator=(const AliAODTracklets& evt);

  virtual ~AliAODTracklets();

  void CreateContainer(Int_t nTracks);
  void DeleteContainer();
  virtual void Clear(Option_t* )         {AliVMultiplicity::Clear(); DeleteContainer();}

  Bool_t SetTracklet(Int_t pos, Double32_t theta, Double32_t phi, Double32_t deltaPhi, Int_t labelL1, Int_t labelL2);


  virtual Int_t    GetNumberOfTracklets() const { return fNTracks; }
  virtual Double_t GetTheta(Int_t i)      const;
  virtual Double_t GetPhi(Int_t i)        const;
  virtual Double_t GetDeltaPhi(Int_t i)   const;
  virtual Int_t    GetLabel(Int_t i, Int_t layer) const;
  virtual void     SetLabel(Int_t i, Int_t layer,Int_t label);
  //
  virtual Double_t* GetTheta()       const {return (Double_t*)fTheta;}
  virtual Double_t* GetPhi()         const {return (Double_t*)fPhi;}
  virtual Double_t* GetDeltPhi()     const {return (Double_t*)fDeltaPhi;}
  virtual Int_t*    GetLabels()      const {return (Int_t*)fLabels;}  
  virtual Int_t*    GetLabels2()     const {return (Int_t*)fLabelsL2;}
  //
  virtual Short_t   GetNumberOfFiredChips(Int_t layer)             const {return fFiredChips[layer];}
  virtual void      SetFiredChips(Int_t layer, Short_t firedChips)       {fFiredChips[layer] = firedChips;}
  //
  virtual void      SetFastOrFiredChips(UInt_t chipKey)                  {fFastOrFiredChips.SetBitNumber(chipKey);}
  virtual const TBits& GetFastOrFiredChips()                       const {return fFastOrFiredChips;}
  virtual void      SetFastOrFiredChipMap(const TBits& fochips)          {fFastOrFiredChips = fochips;}
  virtual Bool_t    TestFastOrFiredChips(UInt_t chipKey)           const {return fFastOrFiredChips.TestBitNumber(chipKey);}
  //
  virtual void   SetFiredChipMap(const TBits & firedChips)               {fClusterFiredChips = firedChips;}
  virtual void   SetFiredChipMap(UInt_t chipKey)                         {fClusterFiredChips.SetBitNumber(chipKey);}
  virtual const  TBits& GetFiredChipMap()                          const {return fClusterFiredChips;}
  virtual Bool_t TestFiredChipMap(UInt_t chipKey)                  const {return fClusterFiredChips.TestBitNumber(chipKey);}

  //
  virtual void Print(Option_t *opt="") const;

  virtual UInt_t GetNumberOfITSClusters(Int_t layer) const { return fITSClusters[layer];}
  virtual void   SetITSClusters(Int_t layer, UInt_t clusters) { fITSClusters[layer] = clusters; }

 protected:
  Int_t      fNTracks;       // Number of tracklets
  Double32_t *fTheta;        //[fNTracks] array with theta values
  Double32_t *fPhi;          //[fNTracks] array with phi values
  Double32_t *fDeltaPhi;     //[fNTracks] array with delta phi values
  Int_t      *fLabels;       //[fNTracks] array with labels of cluster in L1 used for the tracklet
  Int_t      *fLabelsL2;     //[fNTracks] array with labels of cluster in L2 used for the tracklet

  Short_t fFiredChips[2];    // Number of fired chips in the two SPD layers
  UInt_t fITSClusters[6];    // Number of ITS cluster per layer
  TBits fFastOrFiredChips;   // Map of FastOr fired chips
  TBits fClusterFiredChips;  // Map of fired chips (= at least one cluster)

  ClassDef(AliAODTracklets, 5);
};


#endif
