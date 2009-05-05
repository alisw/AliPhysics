#ifndef ALITRDLTUTRACKLET_H
#define ALITRDLTUTRACKLET_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDltuTracklet.h 14830 2006-08-11 17:58:05Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD LTU tracklet                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliTRDltuTracklet : public TObject {
  
 public:

  enum { kNplan = 6 };

  AliTRDltuTracklet();
  AliTRDltuTracklet(Int_t det, Int_t row, Float_t rowz, Float_t slope, Float_t offset
		   ,Float_t time, Int_t ncl, Int_t label, Float_t q);
  virtual         ~AliTRDltuTracklet();

          Bool_t   IsSortable() const { return kTRUE; }
  virtual Int_t    Compare(const TObject *o) const;

          Int_t    GetDetector() const             { return fDetector;                };
          Int_t    GetPlane(Int_t det) const       { return ((Int_t) (det % kNplan)); };
          Int_t    GetRow() const                  { return fRow;                     };
          Int_t    GetNclusters() const            { return fNclusters;               };
          Float_t  GetSlope() const                { return fSlope;                   };
          Float_t  GetOffset() const               { return fY;                       }; 
          Float_t  GetTime0() const                { return fX;                       };
          Float_t  GetRowz() const                 { return fRowz;                    };
          Float_t  GetYproj(Float_t xpl) const;
          Float_t  GetZproj(Float_t xpl) const;
          Int_t    GetLabel() const                { return fLabel;                   };
          Float_t  GetPt(Float_t field) const;
          Float_t  GetQ() const                    { return fQ;                       };

 protected:

          Float_t  fX;                              // Distance vertex to entrance window
          Float_t  fY;                              // Tracklet offset at entrance window
          Float_t  fSlope;                          // Tracklet slope
          Float_t  fRowz;                           // z coordinate of the pad row center
          Int_t    fDetector;                       // Detector number
          Int_t    fRow;                            // Pad row number 
          Int_t    fNclusters;                      // Number of clusters
          Int_t    fLabel;                          // MC track label
          Float_t  fQ;                              // Charge sum divided by number of clusters

  ClassDef(AliTRDltuTracklet,2)                     // LTU tracklet

};

#endif
