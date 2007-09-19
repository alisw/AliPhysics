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

  AliTRDcluster();
  AliTRDcluster(Int_t det, Float_t q, Float_t *pos, Float_t *sig
              , Int_t *tracks, Char_t npads, Short_t *signals
              , UChar_t col, UChar_t row, UChar_t time
              , Char_t timebin, Float_t center, UShort_t volid);
  AliTRDcluster(const AliTRDcluster &c);

  virtual void     AddTrackIndex(Int_t *i); 

          Int_t    IsUsed() const               { return (fQ < 0) ? 1 : 0; }
          void     Use(Int_t = 0)               { fQ = -fQ;                }
    
          Int_t    GetDetector() const          { return fDetector;        }
          Int_t    GetLocalTimeBin() const      { return fLocalTimeBin;    }
          Float_t  GetQ() const                 { return fQ;               }
          Int_t    GetNPads() const             { return fNPads;           }
          Float_t  GetCenter() const            { return fCenter;          }
	  Int_t    GetPadCol() const            { return fPadCol;          }
	  Int_t    GetPadRow() const            { return fPadRow;          }
          Int_t    GetPadTime() const           { return fPadTime;         }
          Short_t *GetSignals()                 { return fSignals;         }
          Float_t  GetSumS() const;

          void     SetLocalTimeBin(Char_t t)    { fLocalTimeBin = t;       }

 protected:
  
          Int_t   fDetector;       //  TRD detector number
          Char_t  fLocalTimeBin;   //  T0-calibrated time bin number
          Float_t fQ;              //  Amplitude 
          Char_t  fNPads;          //  Number of pads in cluster
          Float_t fCenter;         //  Center of the cluster relative to the pad 
	  UChar_t fPadCol;         //  Central pad number in column direction
	  UChar_t fPadRow;         //  Central pad number in row direction
	  UChar_t fPadTime;        //  Uncalibrated time bin number
          Short_t fSignals[7];     //  Signals in the cluster
  
  ClassDef(AliTRDcluster,5)        //  Cluster for the TRD
 
};
#endif
