/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSlat.h,v 1.2 2005/09/19 19:01:09 ivana Exp $

/// \ingroup slat
/// \class AliMpSlat
/// \brief A slat (building block of stations 3, 4 and 5)
/// 
/// Author: Laurent Aphecetche

#ifndef ALI_MP_SLAT_H
#define ALI_MP_SLAT_H

#include <TObject.h>

#ifndef ROOT_TString
#include "TString.h"
#endif

#ifndef ALI_MP_V_SEGMENTATION_H
#include "AliMpVSegmentation.h"
#endif

#ifndef ALI_MP_PAD_H
#include "AliMpPad.h"
#endif

#include "AliMpContainers.h"

class TArrayI;

#ifdef WITH_ROOT
#  include "TExMap.h"
#  include "TObjArray.h"
#else
#  include <vector>
#  include <map>
#endif

class AliMpMotifPosition;
class AliMpPCB;

class AliMpSlat : public TObject
{
 public:

#ifdef WITH_ROOT
  typedef Int_t Size_t;
#else
  typedef UInt_t Size_t;
#endif  
  
  AliMpSlat();
  AliMpSlat(const char* id);
  virtual ~AliMpSlat();

  TVector2 Dimensions() const;

  const char* GetID() const;

  void Add(AliMpPCB* pcbType, const TArrayI& manuList);

  Double_t DX() const;
  Double_t DY() const;

  /// Find the PCB containing the pad at location (ix,any iy).
  AliMpPCB* FindPCB(Int_t ix) const;

  /** Find the index of the PCB containing the pad at location ix.
   Should not be needed except to comply with Sector(), Dpx(), Dpy()
	 interface of old AliMUONVGeometrySegmentation.
	 FIXME: Remove me when VGeometrySegmentation dies at last.
	 */
	Int_t FindPCBIndex(Int_t ix) const;
	
  /// Find the PCB containing location (x,y).
  AliMpPCB* FindPCB(Double_t x, Double_t y) const;

	/** Find the index of the PCB containing the pad at location (x,y).
   Should not be needed except to comply with Sector(), Dpx(), Dpy()
	 interface of old AliMUONVGeometrySegmentation.
	 FIXME: Remove me when VGeometrySegmentation dies at last.
	 */
	Int_t FindPCBIndex(Double_t x, Double_t y) const;

  /// Returns the i-th PCB of this slat.
  AliMpPCB* GetPCB(Size_t i) const;

  /// Returns the MotifPosition containing location (x,y).
  AliMpMotifPosition* FindMotifPosition(Double_t x, Double_t y) const;

  /// Returns the MotifPosition which id is manuid.
  AliMpMotifPosition* FindMotifPosition(Int_t manuid) const;

  /// Returns the MotifPosition containing the pad located at (ix,iy).
  AliMpMotifPosition* FindMotifPosition(Int_t ix, Int_t iy) const;

  /// Returns the number of PCBs of this slat.
  Size_t GetSize() const;

  /// Returns the number of pads in the x-direction contained in this slat.
  Int_t GetNofPadsX() const;
 
  /** Returns the max. number of pads in the x-direction contained in this slat.
      This is a max only as for e.g. non-bending slats, the y-dimension depends
      on the x-position.
  */
  Int_t GetMaxNofPadsY() const;

  void Print(Option_t* option="") const;

 private:
  TString fId;
  Double_t fDX;
  Double_t fDY;
  Int_t fNofPadsX;
  Int_t fMaxNofPadsY;
#ifdef WITH_ROOT
  TObjArray fPCBs; // array of AliMpPCB*
  mutable TExMap fManuMap; // map of int to AliMpMotifPosition*
#else  
  std::vector<AliMpPCB*> fPCBs;
  std::map<int,AliMpMotifPosition*> fManuMap;
#endif
  
  ClassDef(AliMpSlat,1) // A slat for stations 3,4,5
};

#endif
