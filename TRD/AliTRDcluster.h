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

class AliTRDcluster : public AliCluster {

 public:

  enum { kInChamber = BIT(16)};
  enum { kMaskedLeft   = 0
       , kMaskedCenter = 1
       , kMaskedRight  = 2
  };

  AliTRDcluster();
  AliTRDcluster(
    Int_t det, Float_t q, Float_t *pos, Float_t *sig,
    Int_t *tracks, Char_t npads, Short_t *signals,
    UChar_t col, UChar_t row, UChar_t time,
    Char_t timebin, Float_t center, UShort_t volid);
  AliTRDcluster(const AliTRDcluster &c);
  virtual ~AliTRDcluster() {};

  virtual void     AddTrackIndex(Int_t *i); 
  void     Clear(Option_t *o="");
  
  Bool_t   IsEqual(const TObject *otherCluster) const;
  Bool_t   IsInChamber() const             { return TestBit(kInChamber);       }
  Bool_t   IsMasked() const                { return fClusterMasking ? kTRUE : kFALSE; }
  Bool_t   IsShared() const                { return IsClusterShared();}
  Bool_t   IsUsed() const                  { return IsClusterUsed(); }

  UChar_t  GetPadMaskedPosition() const    { return fClusterMasking & 7; }
  UChar_t  GetPadMaskedStatus() const      { return fClusterMasking >> 3; }
  Int_t    GetDetector() const             { return fDetector;      }
  Int_t    GetLocalTimeBin() const         { return fLocalTimeBin;  }
  Float_t  GetQ() const                    { return fQ;             }
  Int_t    GetNPads() const                { return fNPads;         }
  Float_t  GetCenter() const               { return fCenter;        }
  Int_t    GetPadCol() const               { return fPadCol;        }
  Int_t    GetPadRow() const               { return fPadRow;        }
  Int_t    GetPadTime() const              { return fPadTime;       }
  Short_t *GetSignals()                    { return fSignals;       }
  Float_t  GetSumS() const;
  
  void     Print(Option_t* o="") const;

  void     SetLocalTimeBin(Char_t t)       { fLocalTimeBin = t;     }
  void     SetInChamber(Bool_t in = kTRUE) { SetBit(kInChamber,in); }
  void     SetPadMaskedPosition(UChar_t position);
  void     SetPadMaskedStatus(UChar_t status);
  void     SetPadCol(UChar_t inPadCol){ fPadCol = inPadCol;}
  void     SetPadRow(UChar_t inPadRow){ fPadRow = inPadRow;}
  void     SetPadTime(UChar_t inPadTime){ fPadTime = inPadTime;}
  void     SetDetector(Short_t inDetector){ fDetector = inDetector;}
  void     SetQ(Float_t inQ){ fQ = inQ;}
  void     SetClusterMasking(UChar_t inClusterMasking){ fClusterMasking = inClusterMasking;}
  void     SetShared(Bool_t sh  = kTRUE)   { SetBit(AliCluster::kShared,sh);    }
  void     Use(Int_t = 0)                  { SetBit(AliCluster::kUsed, kTRUE);              }

  protected:
    UChar_t fPadCol;         //  Central pad number in column direction
    UChar_t fPadRow;         //  Central pad number in row direction
    UChar_t fPadTime;        //  Uncalibrated time bin number
    Char_t  fLocalTimeBin;   //  T0-calibrated time bin number
    UChar_t fNPads;          //  Number of pads in cluster
    UChar_t fClusterMasking; //  Bit field containing cluster status information;
    Short_t fDetector;       //  TRD detector number
    Short_t fSignals[7];     //  Signals in the cluster
    Float_t fQ;              //  Amplitude 
    Float_t fCenter;         //  Center of the cluster relative to the pad 

    ClassDef(AliTRDcluster, 6)        //  Cluster for the TRD
};
#endif
