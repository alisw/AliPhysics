#ifndef ALIMUONTRIGGERSEGMENTATION_H
#define ALIMUONTRIGGERSEGMENTATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup base
/// \class AliMUONTriggerSegmentation
/// \brief Segmentation for trigger modules
///
/// New implementation of AliMUONVGeometryDESegmentation, based on  
/// mapping package.

#ifndef ALIMUONVGEOMETRYDESEGMENTATION_H
#include "AliMUONVGeometryDESegmentation.h"
#endif

#ifndef ALI_MP_PLANE_TYPE_H
#include "AliMpPlaneType.h"
#endif

#ifndef ALI_MP_PAD_H
#include "AliMpPad.h"
#endif

class AliMpTrigger;
class AliMpTriggerSegmentation;

class AliMUONTriggerSegmentation : public AliMUONVGeometryDESegmentation 
{
 public:

  AliMUONTriggerSegmentation();
  AliMUONTriggerSegmentation(AliMpVSegmentation* segmentation,
                               Int_t detElemId, AliMpPlaneType bendingOrNonBending);
  virtual ~AliMUONTriggerSegmentation();
      
  /// Distance between 1 pad and a position
  virtual Float_t Distance2AndOffset(Int_t /*iX*/, Int_t /*iY*/, 
                                     Float_t /*X*/, Float_t /*Y*/, 
                                     Int_t * /*dummy*/);
  
  virtual Float_t Dpx() const;
  virtual Float_t Dpy() const;
  virtual Float_t Dpx(Int_t sectorCode) const;
  virtual Float_t Dpy(Int_t sectorCode) const;
  
  virtual void Draw(Option_t*/*opt*/ = "");
  
  void FirstPad(Float_t xhit, Float_t yhit, Float_t zhit, 
                Float_t dx, Float_t dy);

  virtual Bool_t HasPad(Float_t /*x*/, Float_t /*y*/, Float_t /*z*/);
  virtual Bool_t HasPad(Int_t ix, Int_t iy);
  
  virtual AliMUONGeometryDirection GetDirection();
  virtual const AliMpVSegmentation* GetMpSegmentation() const;

  virtual Float_t GetAnod(Float_t /*xhit*/) const;

  virtual void GetPadI(Float_t x ,Float_t y ,Int_t   &ix,Int_t &iy);
  virtual void GetPadI(Float_t x, Float_t y , Float_t z, Int_t &ix, Int_t &iy);
  virtual void GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y);
  virtual void GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z);
                  
  virtual void Init(Int_t) {} ///< Not implemented
  
  virtual void IntegrationLimits(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
  
  virtual Int_t ISector();
  virtual Int_t Ix();
  virtual Int_t Iy();
  
  Int_t LineNumber() const;
  
  virtual Int_t MorePads();
 
  virtual void Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, 
                          Int_t Xlist[10], Int_t Ylist[10]);
  virtual void NextPad();
    
  virtual Int_t Npx() const;
  virtual Int_t Npy() const;

  void Print(Option_t* opt="") const;
  
  virtual void SetDAnod(Float_t /*D*/);

  virtual Int_t Sector(Int_t ix, Int_t iy);
  virtual void SetHit(Float_t xhit, Float_t yhit); 
  virtual void SetHit(Float_t xhit, Float_t yhit, Float_t zhit);

  virtual void SetPad(Int_t ix, Int_t iy);
  virtual void SetPadSize(Float_t, Float_t);
  
  virtual void GetNParallelAndOffset(Int_t /*iX*/, Int_t /*iY*/, 
                                     Int_t */*Nparallel*/, Int_t */*Offset*/);
  virtual Int_t SigGenCond(Float_t /*x*/, Float_t /*y*/, Float_t /*z*/);
  virtual void SigGenInit(Float_t /*x*/, Float_t /*y*/, Float_t /*z*/);
  virtual void GiveTestPoints(Int_t &/*n*/, Float_t * /*x*/, Float_t */*y*/) const;
  virtual void SetCorrFunc(Int_t /*dum*/, TF1* /*func*/);
  virtual TF1* CorrFunc(Int_t) const;
  virtual Int_t Sector(Float_t /*x*/, Float_t /*y*/);

public:

    void GetPadLoc2Glo(Int_t ixLoc, Int_t iyLoc, Int_t& ixGlo, Int_t& iyGlo) const;
    void GetPadGlo2Loc(Int_t ixLoc, Int_t iyLoc, Int_t& ixGlo, Int_t& iyGlo) const;
    
    void PC2LA(Int_t ixPC, Int_t iyPC, Int_t& ixLA, Int_t& iyLA) const;
    void LA2PC(Int_t ixLA, Int_t iyLA, Int_t& ixPC, Int_t& iyPC) const;
    
    void IGlo2ILoc(Int_t ixGlo, Int_t iyGlo, Int_t& ixLA, Int_t& iyLA) const;
    void ILoc2IGlo(Int_t ixLA, Int_t iyLA, Int_t& ixGlo, Int_t& iyGlo) const;
    
    Int_t ModuleColNum(Int_t ixGlo) const;
    
protected:

    AliMUONTriggerSegmentation(const AliMUONTriggerSegmentation& rhs);
    AliMUONTriggerSegmentation& operator=(const AliMUONTriggerSegmentation& rhs);

private:
    Int_t fDetElemId;          ///< det elem Id
    AliMpPlaneType fPlaneType; ///< plane type
    const AliMpTrigger* fSlat; ///< slat
    AliMpTriggerSegmentation* fSlatSegmentation; ///< mapping segmentation
//    AliMpVPadIterator* fPadIterator; //!
    AliMpPad fCurrentPad; //!< FIXME: should not be needed, if we externalise the SetPad, SetHit, IntegrationLimits methods which have nothing to do here anyway, together with the iteration methods FirstPad, NextPad, MorePads, which have nothing to do here either.
    Float_t fXhit;        //!< x-position of hit
    Float_t fYhit;        //!< y-position of hit
    Int_t fLineNumber;    ///< Line number of that detection element (from 1 to 9)
    
    ClassDef(AliMUONTriggerSegmentation,1) // Trigger segmentation
};
#endif






