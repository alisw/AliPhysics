/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$ 
// $MpId: AliMpPCB.h,v 1.9 2006/05/24 13:58:24 ivana Exp $ 

/// \ingroup slat
/// \class AliMpPCB
/// \brief A PCB for station 3,4 or 5
/// 
//  Author: Laurent Aphecetche

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

#ifndef ALI_MP_AREA_H
#  include "AliMpArea.h"
#endif

class AliMpSlatMotifMap;
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
      See full doc for the meaning of enveloppe parameters.
  */
  AliMpPCB(AliMpSlatMotifMap* motifMap,
           const char* id, Double_t padSizeX, Double_t padSizeY,
           Double_t enveloppeSizeX, Double_t enveloppeSizeY);
  
  AliMpPCB(const char* id, AliMpMotifSpecial* ms);
  
  AliMpPCB(const AliMpPCB& o);
  AliMpPCB& operator=(const AliMpPCB& o);

  virtual ~AliMpPCB();

  TObject* Clone(const char* newname="") const;

  /** Duplicate this PCB. The copy has the manuids of its motifs changed 
      according to the manuid vector, and its x-offset according to ix 
      and x.
  */    
  AliMpPCB* Clone(const TArrayI& manuids,
                  Int_t ixOffset, Double_t xOffset) const;

  void Copy(TObject& o) const;

  /** Add a motif to this PCB. (ix,iy) are the coordinates of one corner 
    of the motif, in pad-units. Which corner depends on the sign(s) of (ix,iy):
    (ix>0,iy>0) : bottom-left corner
    (ix<0,iy>0) : bottom-right corner
    (ix<0,iy<0) : top-right corner
    (ix>0,iy<0) : top-left corner.
    */
  void Add(AliMpMotifType* motifType, Int_t ix, Int_t iy);

  AliMpArea Area() const;
  
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
  
  Int_t NofPads() const { return fNofPads; }
  
  AliMpSlatMotifMap* MotifMap() const { return fMotifMap; }
  
  void Save() const;
  
 private:
  TString fId; ///< PCB name
  Double_t fPadSizeX; ///< x-size of this PCB's pads (cm)
  Double_t fPadSizeY; ///< y-size of this PCB's pads (cm)
  Double_t fEnveloppeSizeX; ///< max x-size of this PCB (cm)
  Double_t fEnveloppeSizeY; // max y-size of this PCB (cm)
  Double_t fXoffset; ///< x-offset = x of first pad of this PCB (cm)
  Double_t fActiveXmin; ///< min x of an actual pad in this PCB (cm)
  Double_t fActiveXmax; ///< max x of an actual pad in this PCB (cm)
  Int_t fIxmin; ///< min pad index in x
  Int_t fIxmax; ///< max pad index in x
  Int_t fIymin; ///< min pad index in y
  Int_t fIymax; ///< max pad index in y
#ifdef WITH_ROOT
  TObjArray fMotifPositions; ///< array of motifs
#else  
  std::vector<AliMpMotifPosition*> fMotifPositions; ///< array of motif positions
#endif
  Int_t fNofPads; ///< number of pads in this PCB
  AliMpSlatMotifMap* fMotifMap; ///< to keep track of things to avoid duplications of motif and motiftypes, and get proper deletion
  
  ClassDef(AliMpPCB,3) // A PCB for Stations 3,4,5
};

#endif 
