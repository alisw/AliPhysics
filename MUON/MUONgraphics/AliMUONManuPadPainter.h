#ifndef ALIMUONMANUPADPAINTER_H
#define ALIMUONMANUPADPAINTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONManuPadPainter
/// \brief Painter for the pads of one manu
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVPAINTER_H
#  include "AliMUONVPainter.h"
#endif
#ifndef ALI_MP_PAD_H
#  include "AliMpPad.h"
#endif

class AliMUONManuPadPainter : public AliMUONVPainter
{
public:
  AliMUONManuPadPainter(TRootIOCtor* ioCtor);
  AliMUONManuPadPainter();
  AliMUONManuPadPainter(const AliMUONVPainter& mother,
                        Int_t detElemId,
                        Int_t manuId);
  virtual ~AliMUONManuPadPainter();

  /// Clone ourselves
  virtual TObject* Clone(const char* = "") const { return new AliMUONManuPadPainter(*this); }
  
  virtual void ComputeDataRange(const AliMUONVTrackerData& data,
                                Int_t dataIndex,
                                Double_t& dataMin, Double_t& dataMax) const;
    
  virtual char* GetObjectInfo(Int_t px, Int_t py) const;
  
  /// We advertise that we do handle mouse movement
  virtual Bool_t HandleMouseMotion() const { return kTRUE; }
  
  TString NameAtPosition(Double_t x, Double_t y) const;

    virtual TString Describe(const AliMUONVTrackerData& data, Int_t dataIndex,
                             Double_t x=FLT_MAX, Double_t y=FLT_MAX);

  using AliMUONVPainter::PaintArea;
  
  void PaintArea(const AliMUONVTrackerData& data,
                   Int_t dataIndex,
                   Double_t min,
                   Double_t max);
  
  void PaintOutline(Int_t color=-1, Int_t width=-1, Double_t x=FLT_MAX, Double_t y=FLT_MAX);

  /// Whether this painter can be detached from the current view.
  Bool_t CanBeDetached() const { return kFALSE; }

  virtual void DrawHistogramClone(Double_t* values=0x0) const;

  virtual Bool_t IsIncluded() const;
  
private:
    
  void BackupStyle();
  void RestoreStyle();
  AliMpPad PadByPosition(Double_t x, Double_t y) const;
  void PaintPad(const AliMpPad& pad) const;
  
private:
  Int_t fDetElemId; ///< our detection element id
  Int_t fManuId; ///< our manu id
  Int_t fLineColorBck; ///< line color for backup
  Int_t fLineWidthBck; ///< line width for backup
  Int_t fFillColorBck; ///< fill color for backup
  Int_t fFillStyleBck; ///< fill style for backup
  
  ClassDef(AliMUONManuPadPainter,1) // Painter for the pads of one manu
};

#endif
