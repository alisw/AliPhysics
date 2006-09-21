#ifndef ALITRDCLUSTER_H
#define ALITRDCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD cluster                                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliCluster.h"  

class AliTRDrecPoint;

class AliTRDcluster : public AliCluster {

 public:

  AliTRDcluster();
  AliTRDcluster(const AliTRDcluster &c);
  AliTRDcluster(const AliTRDrecPoint &p);

  virtual void     AddTrackIndex(Int_t *i); 

          Int_t    IsUsed() const               { return (fQ < 0) ? 1 : 0;    }
          void     Use(Int_t = 0)               { fQ = -fQ;                   }
  
          Bool_t   From2pad() const             { return TestBit(k2pad);      }
          Bool_t   From3pad() const             { return TestBit(k3pad);      }
          Bool_t   From4pad() const             { return TestBit(k4pad);      }
          Bool_t   From5pad() const             { return TestBit(k5pad);      }
          Bool_t   FromLarge() const            { return TestBit(kLarge);     }
          Bool_t   Isolated() const             { return (TestBit(k2pad) || TestBit(k3pad)); }
 
  virtual void     SetDetector(Int_t d)         { fDetector = d;    } 
          void     SetLocalTimeBin(Int_t t)     { fTimeBin  = t;    }
          void     SetQ(Float_t q)              { fQ        = q;    }
  virtual void     SetX(Float_t x)              { fX        = x;    }
          void     SetCenter(Float_t c)         { fCenter   = c;    }
          void     SetSignals(Short_t *signals);
          void     Set2pad()                    { SetBit(k2pad);  fNPads = 2; }
          void     Set3pad()                    { SetBit(k3pad);  fNPads = 3; }
          void     Set4pad()                    { SetBit(k4pad);  fNPads = 4; }
          void     Set5pad()                    { SetBit(k5pad);  fNPads = 5; }
          void     SetLarge()                   { SetBit(kLarge); fNPads = 6; }
  
  virtual Int_t    GetDetector() const          { return fDetector; }
          Int_t    GetLocalTimeBin() const      { return fTimeBin;  }
          Float_t  GetQ() const                 { return fQ;        }
  virtual Float_t  GetX() const                 { return fX;        }
          Int_t    GetNPads() const             { return fNPads;    }
          Float_t  GetCenter() const            { return fCenter;   }
          Short_t *GetSignals()                 { return fSignals;  }
          Float_t  GetSumS() const;

 protected:

          enum {
            k2pad  = 0x00000001,   //  2 pad cluster
            k3pad  = 0x00000002,   //  3 pad cluster
            k4pad  = 0x00000004,   //  4 pad cluster
            k5pad  = 0x00000008,   //  5 pad cluster
            kLarge = 0x00000016    //  Large cluster
          };
  
          Int_t   fDetector;       //  TRD detector number
          Float_t fX;              //  Local x position (first order proportional to time bin)
          Char_t  fTimeBin;        //  Time bin number within the detector
          Float_t fQ;              //  Amplitude 
          Char_t  fNPads;          //  Number of pads in cluster
          Float_t fCenter;         //  Center of the cluster relative to the pad 
          Short_t fSignals[7];     //  Signals in the cluster
  
  ClassDef(AliTRDcluster,3)        //  Cluster for the TRD
 
};


#endif
