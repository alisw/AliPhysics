#ifndef ALIMUONDEPAINTER_H
#define ALIMUONDEPAINTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONDEPainter
/// \brief A painter for one detection element
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVPAINTER_H
#  include "AliMUONVPainter.h"
#endif

class AliMUONDEPainter : public AliMUONVPainter
{
public:
  AliMUONDEPainter();
  AliMUONDEPainter(TRootIOCtor* ioCtor);
  AliMUONDEPainter(const AliMUONAttPainter& att, Int_t detElemId);
  AliMUONDEPainter(const AliMUONDEPainter& rhs);
  AliMUONDEPainter& operator=(const AliMUONDEPainter& rhs);
  virtual ~AliMUONDEPainter();
  
  /// Clone this object
  virtual TObject* Clone(const char* = "") const { return new AliMUONDEPainter(*this); }

  void ComputeDataRange(const AliMUONVTrackerData& data, Int_t dataIndex, 
                        Double_t& dataMin, Double_t& dataMax) const;
    
  void Copy(TObject& object) const;
  
  /// Return the ID of this detection element
  Int_t DetElemId() const { return fDetElemId; }
  
  void PaintArea(const AliMUONVTrackerData& data, Int_t dataIndex,
                 Double_t min, Double_t max);
    
  TString Describe(const AliMUONVTrackerData& data, Int_t dataIndex,
                   Double_t, Double_t);
    
  virtual AliMUONAttPainter Validate(const AliMUONAttPainter& attributes) const;

  virtual void FillManuList(TObjArray& manuList) const;
  
  virtual Bool_t IsIncluded() const;
  
private:
  Int_t fDetElemId; ///< our id

  ClassDef(AliMUONDEPainter,1) // Detection element painter
};

#endif
