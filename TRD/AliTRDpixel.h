#ifndef ALITRDPIXEL_H
#define ALITRDPIXEL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>

//////////////////////////////////////////////////////
//  Stores the information for one detector pixel   //
//////////////////////////////////////////////////////

class AliTRDpixel : public TObject {

public:

  AliTRDpixel();
  virtual ~AliTRDpixel();

  virtual void    Copy(TObject &p);

  static  Int_t   NTrackPixel()                  { return fgkNTrackPixel; };

  virtual void    SetSignal(Float_t signal)      { fSignal   = signal; };
  virtual void    SetTrack(Int_t i, Int_t track) { fTrack[i] = track;  };

  virtual Float_t GetSignal() const              { return fSignal;     };
  virtual Int_t   GetTrack(Int_t i) const        { return fTrack[i];   };

protected:

  enum { kNTrackPixel = 3 };
  static const Int_t fgkNTrackPixel;       // Maximal number of stored tracks

  Float_t      fSignal;                    // Signal sum
  Int_t        fTrack[kNTrackPixel];       // Tracks contributing to this pixel

  ClassDef(AliTRDpixel,1)                  // Information for one detector pixel   

};

#endif
