#ifndef ALIMUONPCBPAINTER_H
#define ALIMUONPCBPAINTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPCBPainter
/// \brief Implementation of AliMUONVPainter for slat's PCBs
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVPAINTER_H
#  include "AliMUONVPainter.h"
#endif

class AliMUONPCBPainter : public AliMUONVPainter
{
public:
  AliMUONPCBPainter(TRootIOCtor* ioCtor);
  AliMUONPCBPainter();
  AliMUONPCBPainter(const AliMUONAttPainter& att, 
                    Int_t detElemId, 
                    Int_t pcbNumber);
  AliMUONPCBPainter(const AliMUONPCBPainter& rhs);
  AliMUONPCBPainter& operator=(const AliMUONPCBPainter& rhs);

  virtual ~AliMUONPCBPainter();

  /// Clone this object
  virtual TObject* Clone(const char* = "" ) const { return new AliMUONPCBPainter(*this); }
  
  virtual void Copy(TObject& object) const;

  void ComputeDataRange(const AliMUONVTrackerData& data, Int_t dataIndex, 
                        Double_t& dataMin, Double_t& dataMax) const;
  
  TString Describe(const AliMUONVTrackerData& data, Int_t dataIndex,
                   Double_t, Double_t);
    
  using AliMUONVPainter::PaintArea;
  
  void PaintArea(const AliMUONVTrackerData& data, Int_t dataIndex,
                 Double_t min, Double_t max);
  
  Bool_t IsIncluded() const;
  
private:
  Int_t fDetElemId; ///< Detection element this pcb is in
  Int_t fPCBIndex;  ///< Index of this PCB within the detection element
  
  ClassDef(AliMUONPCBPainter,1) // Implementation of AliMUONVPainter for St345 PCBs
};

#endif
