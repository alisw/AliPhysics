#ifndef ALITRDRECPOINT_H
#define ALITRDRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD reconstructed point                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliRecPoint.h"

class AliTRDrecPoint : public AliRecPoint {

 public:

  AliTRDrecPoint();
  AliTRDrecPoint(const char * opt);
  virtual         ~AliTRDrecPoint();

  virtual void     Print(Option_t* ) const {};
  virtual void     AddDigit(Int_t digit);
  virtual void     AddDigit(AliDigitNew& ) {};

          void     SetEnergy(Float_t amp)          { fAmp      = amp;   }
          void     SetDetector(Int_t det)          { fDetector = det;   }
          void     SetLocalPosition(TVector3 &pos);
          void     SetLocalRow(Float_t r)          { fLocPos.SetX(r);   }
          void     SetLocalCol(Float_t c)          { fLocPos.SetY(c);   }
          void     SetLocalTime(Float_t t)         { fLocPos.SetZ(t);   }

          void     SetLocalTimeBin(Int_t tb)       { fTimeBin  = tb;    }
          void     SetTrackingYZ(Float_t fSigmaY = 0.0, Float_t fSigmaZ = 0.0);  

          Int_t    GetDetector() const             { return fDetector;  }
          Int_t    GetDigit(Int_t i = 0) const     { if (i < fMulDigit) return fDigitsList[i]; 
	                                             else               return -1;             }
          Float_t  GetLocalRow() const             { return fLocPos(0); }
          Float_t  GetLocalCol() const             { return fLocPos(1); }
          Float_t  GetLocalTime() const            { return fLocPos(2); }

          Int_t    GetLocalTimeBin() const         { return Int_t(fLocPos(2)); }
          Float_t  GetSigmaY2() const              { return fSigmaY2;   }
          Float_t  GetSigmaZ2() const              { return fSigmaZ2;   }
          Float_t  GetY() const                    { return fY;         }
          Float_t  GetZ() const                    { return fZ;         }

          Int_t    IsUsed() const                  { return fUsed;      }
          void     Use()                           { fUsed++;           }
          Int_t    GetTrackIndex(Int_t i) const    { return fTracks[i]; }
          void     AddTrackIndex(Int_t *i);  

 protected:

          Int_t    fDetector;        //  TRD detector number
          Int_t    fTimeBin;         //  Time bin number within the detector
          Int_t    fUsed;            //  0 initially and incremented if the point is "used"
          Int_t    fTracks[3];       //  Labels of overlapped tracks
          Float_t  fY;               //  Local Rphi coordinate (cm) within tracking sector
          Float_t  fZ;               //  Local Z coordinate (cm) within tracking sector
          Float_t  fSigmaY2;         //  Y variance (cm)
          Float_t  fSigmaZ2;         //  Z variance (cm)  

  ClassDef(AliTRDrecPoint,1)         //  Reconstructed point for the TRD

};

#endif
