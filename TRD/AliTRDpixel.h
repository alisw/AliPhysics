#ifndef TRDpixel_h
#define TRDpixel_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>

//////////////////////////////////////////////////////
//  Stores the information for one detector pixel   //
//////////////////////////////////////////////////////

const Int_t kTrackPixel = 3;

class AliTRDpixel : public TObject {

protected:
  Float_t      fSignal;                    // Signal sum
  Int_t        fTrack[kTrackPixel];        // Tracks contributing to this pixel

public:
  AliTRDpixel();
  virtual ~AliTRDpixel() {};

  virtual void    SetSignal(Float_t signal)      { fSignal   = signal; };
  virtual void    SetTrack(Int_t i, Int_t track) { fTrack[i] = track;  };

  virtual Float_t GetSignal()                    { return fSignal;     };
  virtual Int_t   GetTrack(Int_t i)              { return fTrack[i];   };

  ClassDef(AliTRDpixel,1)

};

#endif
