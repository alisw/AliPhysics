/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$ 
// $MpId: AliMpPCB.h,v 1.4 2005/09/19 19:01:09 ivana Exp $ 

/// \ingroup slat
/// \class AliMpPCB
/// \brief A PCB for station 3,4 or 5
/// 
/// Author: Laurent Aphecetche

#ifndef ALIMPPCB_H
#define ALIMPPCB_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ROOT_TString
#  include "TString.h"
#endif

#  ifndef ROOT_TArraI
#    include "TArrayI.h"
#  endif

#include "AliMpContainers.h"

#ifdef WITH_ROOT
#  ifndef ROOT_TObjArray
#    include "TObjArray.h"
#  endif
#else
#  include <vector>
#endif

class AliMpMotifPosition;
class AliMpMotifType;
class AliMpMotifSpecial;

class AliMpPCB : public TObject
{
 public:

#ifdef WITH_ROOT
  typedef Int_t  Size_t;
#else  
  typedef UInt_t Size_t;
#endif
  
  AliMpPCB();
  /** Ctor. The sizes are given in mm.
      enveloppe is due to the fact that not all PCBs are "full" of pads,
      e.g. the rounded or short ones miss some pads, but the enveloppe 
      is a virtual size that should be constant across the slats, 
      and is 400x400 mm.
      It's a usefull notion to compute e.g. slat center in a uniform way, 
      considering that a slat is N PCBs, of the same "virtual" size, that of 
      the enveloppe.
  */
  AliMpPCB(const char* id, Double_t padSizeX, Double_t padSizeY,
	   Double_t enveloppeSizeX, Double_t enveloppeSizeY);
  
  AliMpPCB(const char* id, AliMpMotifSpecial* ms);
  
  AliMpPCB(const AliMpPCB&);
  AliMpPCB& operator=(const AliMpPCB&);

  virtual ~AliMpPCB();

  TObject* Clone(const char* newname="") const;

  /** Duplicate this PCB. The copy has the manuids of its motifs changed 
      according to the manuid vector, and its x-offset according to ix 
      and x.
  */    
  AliMpPCB* Clone(const TArrayI& manuids,
                  Int_t ixOffset, Double_t xOffset) const;

  void Copy(TObject&) const;

  /** Add a motif to this PCB. (ix,iy) are the coordinates of one corner 
    of the motif, in pad-units. Which corner depends on the sign(s) of (ix,iy):
    (ix>0,iy>0) : bottom-left corner
    (ix<0,iy>0) : bottom-right corner
    (ix<0,iy<0) : top-right corner
    (ix>0,iy<0) : top-left corner.
    */
  void Add(AliMpMotifType* motifType, Int_t ix, Int_t iy);

  void Print(Option_t* option = "") const;

  Double_t ActiveDX() const;
  Double_t ActiveDY() const;

  Double_t DX() const;
  Double_t DY() const;

  Double_t X() const;
  Double_t Y() const;

  Double_t Xmin() const;
  Double_t Xmax() const;
 
  Double_t ActiveXmin() const;
  Double_t ActiveXmax() const;

  Double_t Ymin() const;
  Double_t Ymax() const;

  Double_t PadSizeX() const;
  Double_t PadSizeY() const;

  /** Returns the i-th motifPosition of this PCB.
      i : [0..GetSize()-1]
  */
  AliMpMotifPosition* GetMotifPosition(Size_t i) const;

  /// Returns the motifPosition which contains the pad at (ix,iy).
  AliMpMotifPosition* FindMotifPosition(Int_t ix, Int_t iy) const;

  /// Returns the motifPosition which contains the pad at (x,y).
  AliMpMotifPosition* FindMotifPosition(Double_t x, Double_t y) const;

  /// The number of motifs, aka manus.
  Size_t GetSize() const;

  Int_t GetNofPadsX() const;
  Int_t GetNofPadsY() const;

  Int_t Ixmin() const;
  Int_t Ixmax() const;
  
  Int_t Iymin() const;
  Int_t Iymax() const;
  
  const char* GetID() const;
  
 private:
  TString fId;
  Double_t fPadSizeX;
  Double_t fPadSizeY;
  Double_t fEnveloppeSizeX;
  Double_t fEnveloppeSizeY;
  Double_t fXoffset;
  Double_t fActiveXmin;
  Double_t fActiveXmax;
  Int_t fIxmin;
  Int_t fIxmax;
  Int_t fIymin;
  Int_t fIymax;
#ifdef WITH_ROOT
  TObjArray fMotifs;
#else  
  std::vector<AliMpMotifPosition*> fMotifs;
#endif

  ClassDef(AliMpPCB,1) // A PCB for Stations 3,4,5
};

#endif 
