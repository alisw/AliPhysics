#ifndef ALITRDCLUSTER_H
#define ALITRDCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
 
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD cluster                                                              //
//                                                                           //
/////////////////////////////////////////////////////////////////////////////// 

class AliTRDrecPoint;

class AliTRDcluster : public TObject {

 public:

  AliTRDcluster();
  AliTRDcluster(const AliTRDcluster &c);
  AliTRDcluster(const AliTRDrecPoint &p);
  virtual ~AliTRDcluster();
  AliTRDcluster &operator=(const AliTRDcluster &c);

  virtual void    Copy(TObject &c);
  virtual void    AddTrackIndex(Int_t *i);
   
          Int_t   GetDetector() const             { return fDetector; }
          Int_t   GetLocalTimeBin() const         { return fTimeBin;  }

          Float_t GetSigmaY2() const              { return fSigmaY2; }
          Float_t GetSigmaZ2() const              { return fSigmaZ2; }
          Float_t GetY() const                    { return fY; }
          Float_t GetZ() const                    { return fZ; }
          Float_t GetQ() const                    { return fQ; }

          Int_t   IsUsed() const                  { return (fQ < 0) ? 1 : 0; }
          void    Use()                           { fQ = -fQ; }
          Int_t   GetTrackIndex(Int_t i) const    { return fTracks[i]; }

          Bool_t  From2pad() const                { return TestBit(k2pad);  }
          Bool_t  From3pad() const                { return TestBit(k3pad);  }
          Bool_t  From4pad() const                { return TestBit(k4pad);  }
          Bool_t  From5pad() const                { return TestBit(k5pad);  }
          Bool_t  FromLarge() const               { return TestBit(kLarge); }
          Bool_t  Isolated() const                { return (TestBit(k2pad) || TestBit(k3pad)); }
 
          void    SetDetector(Int_t d)            { fDetector  = d; }
          void    SetLocalTimeBin(Int_t t)        { fTimeBin   = t; }
          void    SetQ(Float_t q)                 { fQ         = q; }
          void    SetY(Float_t y)                 { fY         = y; }
          void    SetZ(Float_t z)                 { fZ         = z; }
          void    SetTrackIndex(Int_t i, Int_t t) { fTracks[i] = t; } 
          void    SetSigmaY2(Float_t s)           { fSigmaY2   = s; }
          void    SetSigmaZ2(Float_t s)           { fSigmaZ2   = s; }

          void    Set2pad()                       { SetBit(k2pad);  }
          void    Set3pad()                       { SetBit(k3pad);  }
          void    Set4pad()                       { SetBit(k4pad);  }
          void    Set5pad()                       { SetBit(k5pad);  }
          void    SetLarge()                      { SetBit(kLarge); }

 protected:

  enum {
    k2pad  = 0x00000001,   // 2 pad cluster
    k3pad  = 0x00000002,   // 3 pad cluster
    k4pad  = 0x00000004,   // 4 pad cluster
    k5pad  = 0x00000008,   // 5 pad cluster
    kLarge = 0x00000016    // Large cluster
  };

  Int_t    fDetector;       // TRD detector number
  Int_t    fTimeBin;        // Time bin number within the detector

  Int_t    fTracks[3];      // labels of overlapped tracks
  Float_t  fQ;              // amplitude 
  Float_t  fY;              // local Rphi coordinate (cm) within tracking sector
  Float_t  fZ;              // local Z coordinate (cm) within tracking sector
  Float_t  fSigmaY2;        // Y variance (cm)
  Float_t  fSigmaZ2;        // Z variance (cm)  

  ClassDef(AliTRDcluster,1) // Cluster for the TRD
 
};

#endif
