#ifndef ALIMUONCHAMBERPAINTER_H
#define ALIMUONCHAMBERPAINTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONChamberPainter
/// \brief Painter for one (plane of one) chamber
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVPAINTER_H
#  include "AliMUONVPainter.h"
#endif
#ifndef ALI_MP_CATHOD_TYPE_H
#  include "AliMpCathodType.h"
#endif
#ifndef ALI_MP_PLANE_TYPE_H
#  include "AliMpPlaneType.h"
#endif

class AliMUONChamberPainter : public AliMUONVPainter
{
public:
  AliMUONChamberPainter();
  AliMUONChamberPainter(TRootIOCtor* ioCtor);
  AliMUONChamberPainter(const AliMUONAttPainter& att, Int_t chamberId);
  AliMUONChamberPainter(const AliMUONChamberPainter& rhs);
  AliMUONChamberPainter& operator=(const AliMUONChamberPainter& rhs);
  
  virtual ~AliMUONChamberPainter();

  void ComputeDataRange(const AliMUONVTrackerData& data, Int_t dataIndex, 
                        Double_t& dataMin, Double_t& dataMax) const;
    
  /// Clone ourselves
  virtual TObject* Clone(const char* = "") const { return new AliMUONChamberPainter(*this); }

  virtual void Copy(TObject& object) const;
  
  using AliMUONVPainter::PaintArea;
  
  void PaintArea(const AliMUONVTrackerData& data, Int_t dataIndex,
                 Double_t min, Double_t max);
    
  TString Describe(const AliMUONVTrackerData& data, Int_t dataIndex,
                   Double_t, Double_t);
    
  AliMUONAttPainter Validate(const AliMUONAttPainter& attributes) const;

  Bool_t IsIncluded() const;
  
private:
  Int_t fChamberId; ///< our identifier (0..n)
  
  ClassDef(AliMUONChamberPainter,1) // Painter for one chamber
};

#endif
