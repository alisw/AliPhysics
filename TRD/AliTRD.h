#ifndef TRD_H
#define TRD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set: TRD     //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h" 
#include "AliDigit.h"
#include "AliTRDconst.h"

//_____________________________________________________________________________
class AliTRD : public AliDetector {
 
protected:
  Int_t         fGasMix;                         // Gas mixture. 0: Xe/Isobutane 1: Xe/CO2

  Float_t       fClengthI[kNplan];               // Length of the inner chambers
  Float_t       fClengthM1[kNplan];              // Length of the middle chambers
  Float_t       fClengthM2[kNplan];              // 
  Float_t       fClengthO1[kNplan];              // Length of the outer chambers
  Float_t       fClengthO2[kNplan];              // 
  Float_t       fClengthO3[kNplan];              // 
  Float_t       fCwidth[kNplan];                 // Width of the chambers

  Int_t         fRowMax[kNplan][kNcham][kNsect]; // Number of pad-rows
  Int_t         fColMax[kNplan];                 // Number of pad-columns
  Int_t         fTimeMax;                        // Number of time buckets

  Float_t       fRow0[kNplan][kNcham][kNsect];   // Row-position of pad 0
  Float_t       fCol0[kNplan];                   // Column-position of pad 0
  Float_t       fTime0[kNplan];                  // Time-position of pad 0

  Float_t       fRowPadSize;                     // Pad size in z-direction
  Float_t       fColPadSize;                     // Pad size in rphi-direction
  Float_t       fTimeBinSize;                    // Size of the time buckets

  Int_t         fHole;                           // Geometry without (0) / with (1) hole 

  TClonesArray *fClusters;                       // List of clusters
  Int_t         fNclusters;                      // Number of clusters

public:
  AliTRD();
  AliTRD(const char *name, const char *title);
  virtual        ~AliTRD();
  virtual void    AddHit(Int_t, Int_t*, Float_t*);
  virtual void    AddDigit(Int_t*, Int_t*);    
  virtual void    AddCluster(Int_t*, Int_t*, Float_t*);    
  virtual void    BuildGeometry();
  virtual void    CreateGeometry();
  virtual void    CreateMaterials();
  virtual void    DrawModule();
  Int_t           DistancetoPrimitive(Int_t px, Int_t py);
  TClonesArray   *Cluster()             { return fClusters; };
  virtual void    Init();
  virtual Int_t   IsVersion() const = 0;
  virtual void    MakeBranch(Option_t* option);     
  virtual void    StepManager() = 0; 
  virtual void    SetTreeAddress();

  virtual void    SetHits(Int_t )       {};
  virtual void    SetSensPlane(Int_t)   {};
  virtual void    SetSensChamber(Int_t) {};
  virtual void    SetSensSector(Int_t ) {};

  virtual void    SetGasMix(Int_t imix = 0);

  virtual void    SetRowPadSize(Float_t size)          { fRowPadSize    = size; };
  virtual void    SetColPadSize(Float_t size)          { fColPadSize    = size; };
  virtual void    SetTimeBinSize(Float_t size)         { fTimeBinSize   = size; };

  virtual Float_t GetRowPadSize()                      { return fRowPadSize;  };
  virtual Float_t GetColPadSize()                      { return fColPadSize;  };
  virtual Float_t GetTimeBinSize()                     { return fTimeBinSize; };

  virtual Int_t   GetRowMax(Int_t p, Int_t c, Int_t s) { return fRowMax[p-1][c-1][s-1]; };
  virtual Int_t   GetColMax(Int_t p)                   { return fColMax[p-1];           };
  virtual Int_t   GetTimeMax()                         { return fTimeMax;               };

  virtual Int_t   Hole()                               { return fHole; };

  ClassDef(AliTRD,1)                // Transition Radiation Detector base class

};

//_____________________________________________________________________________
class AliTRDhit : public AliHit {

public:
  Int_t        fSector;     // TRD sector number
  Int_t        fChamber;    // TRD chamber number
  Int_t        fPlane;      // TRD plane number 
  Float_t      fQ;          // Charge created by a hit (slow simulator only)
 
public:
  AliTRDhit() {}
  AliTRDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliTRDhit() {};
 
  ClassDef(AliTRDhit,1)     // Hits for Transition Radiation Detector

};

//_____________________________________________________________________________
class AliTRDdigit : public AliDigit {

public:
  Int_t        fSector;     // TRD sector number
  Int_t        fChamber;    // TRD chamber number
  Int_t        fPlane;      // TRD plane number
  Int_t        fRow;        // Pad row number
  Int_t        fCol;        // Pad col number
  Int_t        fTime;       // Time bucket
  Int_t        fSignal;     // Signal amplitude

public:
  AliTRDdigit() {};
  AliTRDdigit(Int_t *tracks, Int_t *digits);
  virtual ~AliTRDdigit() {};

  ClassDef(AliTRDdigit,1)   // Digits for Transition Radiation Detector

};

//_____________________________________________________________________________
class AliTRDcluster : public TObject {

public:
  Int_t    fSector;         // TRD sector number
  Int_t    fChamber;        // TRD chamber number
  Int_t    fPlane;          // TRD plane number

  Int_t    fTimeSlice;      // Timeslice in chamber where cluster has been found
  Int_t    fEnergy;         // Charge sum of this cluster

  Float_t  fX;              // X coord in ALICE reference frame
  Float_t  fY;              // Y coord in ALICE reference frame
  Float_t  fZ;              // Z coord in ALICE reference frame

  Int_t    fTracks[3];      // Track information

public:
  AliTRDcluster() {};
  AliTRDcluster(Int_t *tracks, Int_t *cluster, Float_t *pos);
  virtual ~AliTRDcluster() {};

  inline virtual Int_t *GetTracks() { return &fTracks[0]; }

  ClassDef(AliTRDcluster,1) // Cluster for Transition Radiation Detector

};

#endif
