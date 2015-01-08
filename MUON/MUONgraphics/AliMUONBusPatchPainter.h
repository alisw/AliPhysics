#ifndef ALIMUONBUSPATCHPAINTER_H
#define ALIMUONBUSPATCHPAINTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONBusPatchPainter
/// \brief A painter for one buspatch
/// 
// Author Laurent Aphecetche, Subatech

#ifndef AliMUONVPAINTER_H
#  include "AliMUONVPainter.h"
#endif
#ifndef ALI_MP_PLANE_TYPE_H
#  include "AliMpPlaneType.h"
#endif

class AliMUONBusPatchPainter : public AliMUONVPainter
{
public:
  
  AliMUONBusPatchPainter();
  AliMUONBusPatchPainter(TRootIOCtor* ioCtor);
  AliMUONBusPatchPainter(const AliMUONAttPainter& att, Int_t busPatchId);
  AliMUONBusPatchPainter(const AliMUONBusPatchPainter& rhs);
  AliMUONBusPatchPainter& operator=(const AliMUONBusPatchPainter& rhs);
  virtual ~AliMUONBusPatchPainter();
  
  /// Clone ourselves
  virtual TObject* Clone(const char* = "") const { return new AliMUONBusPatchPainter(*this); }
  
  void ComputeDataRange(const AliMUONVTrackerData& data, Int_t dataIndex, 
                        Double_t& dataMin, Double_t& dataMax) const;
    
  virtual void Copy(TObject& object) const;
  
  using AliMUONVPainter::PaintArea;
  
  void PaintArea(const AliMUONVTrackerData& data, Int_t dataIndex,
                 Double_t min, Double_t max);
    
  TString Describe(const AliMUONVTrackerData& data, Int_t dataIndex, 
                   Double_t x = FLT_MAX, Double_t y = FLT_MAX);

  virtual AliMUONAttPainter Validate(const AliMUONAttPainter& attributes) const;

  virtual Bool_t IsIncluded() const;
  
private:
  Int_t fBusPatchId; ///< our identifier
  
  ClassDef(AliMUONBusPatchPainter,1) // Painter for one buspatch
};

#endif
