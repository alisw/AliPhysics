#ifndef ALIMUONMANUPAINTER_H
#define ALIMUONMANUPAINTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONManuPainter
/// \brief Painter for one manu (not the pads, only the manu)
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVPAINTER_H
#  include "AliMUONVPainter.h"
#endif

class AliMUONManuPainter : public AliMUONVPainter
{
public:

  AliMUONManuPainter(TRootIOCtor* ioCtor);
  AliMUONManuPainter();
  AliMUONManuPainter(const AliMUONAttPainter& att, 
                     Int_t detElemId, 
                     Int_t manuId);
  AliMUONManuPainter(const AliMUONManuPainter& rhs);
  AliMUONManuPainter& operator=(const AliMUONManuPainter& rhs);
  
  virtual ~AliMUONManuPainter();
  
  virtual void ComputeDataRange(const AliMUONVTrackerData& data, Int_t dataIndex, 
                                Double_t& dataMin, Double_t& dataMax) const;
    
  /// Clone ourselves
  virtual TObject* Clone(const char* = "") const { return new AliMUONManuPainter(*this); }

  virtual void Copy(TObject& object) const;
  
    virtual TString Describe(const AliMUONVTrackerData& data, Int_t dataIndex,
                             Double_t x=FLT_MAX, Double_t y=FLT_MAX);

    void PaintArea(const AliMUONVTrackerData& data, Int_t dataIndex,
                   Double_t min, Double_t max);
  
    virtual AliMUONAttPainter Validate(const AliMUONAttPainter& attributes) const;

    virtual void FillManuList(TObjArray& manuList) const;
    
    virtual Bool_t IsIncluded() const;
    
private:
  Int_t fDetElemId; ///< our detection element id
  Int_t fManuId; ///< our manu id
  
  ClassDef(AliMUONManuPainter,1) // Painter for one manu (not the pads, only the manu)
};

#endif
