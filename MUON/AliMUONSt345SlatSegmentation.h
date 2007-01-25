#ifndef ALIMUONST345SLATSEGMENTATION_H
#define ALIMUONST345SLATSEGMENTATION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup base
/// \class AliMUONSt345SlatSegmentation
/// \brief Segmentation for slat modules

#ifndef ALIMUONVGEOMETRYDESEGMENTATION_H
#include "AliMUONVGeometryDESegmentation.h"
#endif

#ifndef ALI_MP_PLANE_TYPE_H
#include "AliMpPlaneType.h"
#endif

#ifndef ALI_MP_PAD_H
#include "AliMpPad.h"
#endif

class AliMpSlat;
class AliMpSlatSegmentation;
class AliMpVPadIterator;

class AliMUONSt345SlatSegmentation : public AliMUONVGeometryDESegmentation
{
 public:

  AliMUONSt345SlatSegmentation();
  AliMUONSt345SlatSegmentation(AliMpVSegmentation* segmentation,
                                 Int_t detElemId,
				 AliMp::PlaneType bendingOrNonBending);
  virtual ~AliMUONSt345SlatSegmentation();

  void FirstPad(Float_t xhit, Float_t yhit, Float_t zhit, 
		Float_t dx, Float_t dy);

  void  NextPad();

  Int_t MorePads();

  Int_t  Ix();
  Int_t  Iy();
  Int_t  ISector();

  Float_t Distance2AndOffset(Int_t ix, Int_t iy, 
			     Float_t X, Float_t Y, 
			     Int_t* dummy);

  void GetNParallelAndOffset(Int_t iX, Int_t iY,
			     Int_t* Nparallel, Int_t* Offset);

  void Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, 
		  Int_t Xlist[10], Int_t Ylist[10]);

  Int_t  Sector(Int_t ix, Int_t iy);
  Int_t  Sector(Float_t x, Float_t y);

  void  IntegrationLimits(Float_t& x1, Float_t& x2, 
			  Float_t& y1, Float_t& y2);

  void GiveTestPoints(Int_t& n, Float_t* x, Float_t* y) const;

  Int_t SigGenCond(Float_t x, Float_t y, Float_t z);
  void SigGenInit(Float_t x, Float_t y, Float_t z);

  void SetCorrFunc(Int_t isec,  TF1* func);

  TF1* CorrFunc(Int_t isec) const;

  void SetPadSize(float x,float y);

  void SetDAnod(float d);

  void Init(int /*chamber*/) {}  ///< Not implemented
  void Draw(Option_t* opt = "");

  Float_t Dpx() const;
  Float_t Dpy() const;
  
  Float_t Dpx(int ipcb) const;
  Float_t Dpy(int ipcb) const;
  
  Float_t GetAnod(Float_t xhit) const;

  Int_t Npx() const;
  Int_t Npy() const;

  /// Sets the current pad.
  void SetPad(Int_t ix,Int_t iy);
  
  /// Sets the current hit.
  void SetHit(Float_t x, Float_t y, Float_t zIsNotUsed);
 
  AliMUONGeometryDirection GetDirection();// { return kDirUndefined; }

  const AliMpVSegmentation* GetMpSegmentation() const;

  /// \deprecated. Use the one below w/o z instead.
  void GetPadC(Int_t ix, Int_t iy, Float_t& x, Float_t& y, Float_t& z);

  /// From pad indices to coordinates (cm).
  void GetPadC(Int_t ix, Int_t iy, Float_t& x, Float_t& y);

  // Transform from pad to real coordinates
  void GetPadI(Float_t x ,Float_t y ,Int_t   &ix,Int_t &iy);

  // to be deprecated. Use the one above w/o z instead.
  void GetPadI(Float_t x, Float_t y , Float_t z, Int_t &ix, Int_t &iy);


  /// Whether a pad exists at a given position.
  Bool_t HasPad(Float_t x, Float_t y, Float_t z);

  /// Whether a pad exists, given its indices.
  Bool_t HasPad(Int_t ix, Int_t iy);
  
  /// Print.
  void Print(Option_t* opt = "") const;

 protected:
  AliMUONSt345SlatSegmentation(const AliMUONSt345SlatSegmentation& right);
  AliMUONSt345SlatSegmentation&  operator = (const AliMUONSt345SlatSegmentation& right);
     
 private:

  Int_t fDetElemId;                ///< det element Id
	AliMp::PlaneType fPlaneType; ///< plane type
  const AliMpSlat* fSlat;          ///< slat
  AliMpSlatSegmentation* fSlatSegmentation; ///< slat segmentation
  AliMpVPadIterator* fPadIterator; //!< pad iterator
  AliMpPad fCurrentPad; //!< FIXME: should not be needed, if we externalise the SetPad, SetHit, IntegrationLimits methods which have nothing to do here anyway, together with the iteration methods FirstPad, NextPad, MorePads, which have nothing to do here either.
  Float_t fXhit;        //!<  x-position of hit
  Float_t fYhit;        //!<  y-position of hit
  ClassDef(AliMUONSt345SlatSegmentation,4) // St345 segmentation
};

#endif
