#ifndef ALITRDCLUSTERMI_H
#define ALITRDCLUSTERMI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD cluster, alternative version                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDcluster.h"  
#include "TMath.h"  

class AliTRDrecPoint;

class AliTRDclusterMI : public AliTRDcluster {

 public:

  AliTRDclusterMI();
  AliTRDclusterMI(const AliTRDcluster &c);
  AliTRDclusterMI(const AliTRDrecPoint &p);

          void     SetRmsY(Float_t rmsy)          { fRmsY   = rmsy;                   }
          void     SetNPads(Int_t npads)          { fNPads  = npads;                  }
          void     SetRelPos(Float_t pos)         { fRelPos = TMath::Nint(pos*128.0); }

          Float_t  GetRmsY() const                { return fRmsY;                     }
          Char_t   GetNPads() const               { return fNPads;                    }
          Float_t  GetRelPos() const              { return float(fRelPos)/128.0;      }

 protected:

          Float_t  fRmsY;                         // RMS in y direction ????
          Char_t   fNPads;                        // Number of pads ????
          Char_t   fRelPos;		       	  // Relative position ????

  ClassDef(AliTRDclusterMI,2)                     // ClusterMI for the TRD
 
};

#endif
