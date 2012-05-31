#ifndef ALITRDCHAMBERTIMEBIN_H
#define ALITRDCHAMBERTIMEBIN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id: AliTRDchamberTimeBin.h 22646 2007-11-29 18:13:40Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  A TRD layer in a single stack                                         //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class AliTRDcluster;
class AliTRDReconstructor;

class AliTRDchamberTimeBin : public TObject
{
public:
  enum AliTRDchamberTimeBinStatus {
     kT0    = BIT(14) // is the T0 time bin
    ,kOwner = BIT(15) // is owner of the clusters
  };
  enum{
     kMaxClustersLayer = 150
    ,kMaxRows = 16
  };

  AliTRDchamberTimeBin(Int_t plane=-1, Int_t stack=-1, Int_t sector=-1, Double_t z0=-1., Double_t zLength=-1.);
  //AliTRDchamberTimeBin(const AliTRDpropagationLayer &layer, Double_t z0, Double_t zLength, UChar_t stackNr);
  //AliTRDchamberTimeBin(const AliTRDpropagationLayer &layer);
  AliTRDchamberTimeBin(const AliTRDchamberTimeBin &layer);
  ~AliTRDchamberTimeBin();
  operator       Int_t() const  { return fN; }
  AliTRDchamberTimeBin&
                 operator=(const AliTRDchamberTimeBin &myLayer);
  AliTRDcluster* operator[](const Int_t i) const {
    return ((i < fN) && (i >= 0)) ? fClusters[i] : NULL;
  }

  void           Bootstrap(const AliTRDReconstructor *rec, Int_t det);
  void           BuildIndices(Int_t iter = 0);
  void           BuildCond(AliTRDcluster * const cl, Double_t *cond, UChar_t Layer, Double_t theta=0., Double_t phi=0.);
  void           Clear(const Option_t *opt = NULL);
  AliTRDcluster* GetCluster(Int_t index) const {return index < fN && index >= 0 ? fClusters[index] : NULL;}
  Int_t          GetGlobalIndex(Int_t index) const {return ((index < fN) && (index >= 0)) ? fIndex[index] : 0; }
  void           GetClusters(const Double_t * const cond, Int_t *index, Int_t& ncl, Int_t BufferSize = kMaxClustersLayer);
  AliTRDcluster* GetNearestCluster(Double_t *cond);
  Double_t       GetX()                            const {
  return fX;      }
  Double_t       GetZ0()                           const { return fZ0;     }
  Double_t       GetDZ0()                          const { return fZLength;}
  Int_t          GetNClusters()                    const { return fN; }
  Int_t          GetPlane()                        const { return fPlane;  }
  Int_t          GetStack()                        const { return fStack;  }
  Int_t          GetSector()                       const { return fSector; }
  void           InsertCluster(AliTRDcluster *c, UInt_t index);

  Bool_t         IsT0() const {return TestBit(kT0);}
  Bool_t         IsOwner() const {return TestBit(kOwner);}

  void           Print(Option_t *opt=NULL) const;
  Int_t          SearchNearestCluster(Double_t y, Double_t z, Double_t Roady, Double_t Roadz) const;
  void           SetRange(Float_t z0, Float_t zLength);
  void           SetNRows(Int_t nRows){ fNRows = nRows; }
  void           SetPlane(Int_t plane){ fPlane = plane; }
  void           SetReconstructor(const AliTRDReconstructor *rec) {fkReconstructor = rec;}
  void           SetStack(Int_t stack){ fStack = stack; }
  void           SetSector(Int_t sector){ fSector = sector; }
  void           SetOwner(Bool_t copy=kTRUE);
  void           SetT0(Bool_t set=kTRUE) {SetBit(kT0, set);}
  void           SetX(Double_t x) {fX = x;}
private:
  void           Copy(TObject &o) const;
  Int_t          Find(Float_t y) const;
  Int_t          FindYPosition(Double_t y, UChar_t z, Int_t nClusters) const;
  Int_t          FindNearestYCluster(Double_t y, UChar_t z) const;

private:
  const AliTRDReconstructor *fkReconstructor; //! Global TRD reconstructor
  Char_t        fPlane;                       //! Plane number
  Char_t        fStack;                       //! Stack number in supermodule
  Char_t        fSector;                      //! Sector mumber
  Char_t        fNRows;                       //! Number of pad rows in the chamber
  UChar_t       fPositions[kMaxRows];         //! Starting index of clusters in pad row
  Int_t         fN;                           //! Number of clusters
  AliTRDcluster *fClusters[kMaxClustersLayer];//  Array of pointers to clusters
  UInt_t        fIndex[kMaxClustersLayer];    //! Array of cluster indexes
  Double_t      fX;                           //! Radial position of tb

  // obsolete !!
  Double_t      fZ0;                          //  Starting position of the layer in Z direction
  Double_t      fZLength;                     //  Length of the layer in Z direction

  ClassDef(AliTRDchamberTimeBin, 2)           //  Tracking propagation layer for one time bin in chamber
};

#endif	// ALITRDCHAMBERTIMEBIN_H

