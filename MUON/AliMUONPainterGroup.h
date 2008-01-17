#ifndef ALIMUONPAINTERGROUP_H
#define ALIMUONPAINTERGROUP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterGroup
/// \brief A group of AliMUONVPainter
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONVPainter;
class AliMUONVTrackerData;

class AliMUONPainterGroup : public TObject
{
public:
  AliMUONPainterGroup();
  AliMUONPainterGroup(const char* type, Int_t depth);
  virtual ~AliMUONPainterGroup();
  
  Bool_t Add(AliMUONVPainter* painter);
  
  void ComputeDataRange(Double_t& dataMin, Double_t& dataMax);
  
  /// Return the data we are plotting
  AliMUONVTrackerData* Data() const { return fData; }
  
  /// Return the index within the data
  Int_t DataIndex() const { return fDataIndex; }
  
  /// Max data we are plotting
  Double_t DataMax() const { return fDataMax; }
  
  /// Min data we are plotting
  Double_t DataMin() const { return fDataMin; }
  
  /// Depth
  Int_t Depth() const { return fDepth; }
  
  void Draw(Option_t* opt="");
  
  AliMUONVPainter* First() const;
  
  /// We are sortable (by type)
  Bool_t IsSortable() const { return kTRUE; }
  
  Int_t Compare(const TObject* obj) const;
  
  /// Whether we should outline ourselves
  Bool_t IsOutlined() const { return fIsOutlined; }
  
  /// Whether we are the plotting group
  Bool_t IsPlotter() const { return fData != 0 && fDataIndex >= 0; }

  /// Whether we are the responder group
  Bool_t IsResponder() const { return fIsResponder; }
  
  /// Whether we are visible
  Bool_t IsVisible() const { return fIsVisible; }
  
  Bool_t Matches(const char* pattern) const;
  
  void Print(Option_t* opt="") const;
  
  void SetData(AliMUONVTrackerData* data, Int_t dataIndex);
  
  /// Set the data range
  void SetDataRange(Double_t min, Double_t max)
  { fDataMin = min; fDataMax = max; }
  
  Int_t GetLineColor() const;
  
  Int_t GetLineWidth() const;
  
  void SetLine(Int_t lineColor, Int_t lineWidth);
  
  /// Set the outlined flag
  void SetOutlined(Bool_t flag=kTRUE) { fIsOutlined = flag; }
  
  /// Set the responder flag
  void SetResponder(Bool_t flag=kTRUE) { fIsResponder = flag; }
  
  /// Set the visible flag
  void SetVisible(Bool_t flag=kTRUE) { fIsVisible = flag; }
  
  /// Our type
  const char* Type() const { return fType.Data(); }
  
private:
  /// Not implemented
  AliMUONPainterGroup(const AliMUONPainterGroup& rhs);
  /// Not implemented
  AliMUONPainterGroup& operator=(const AliMUONPainterGroup& rhs);
  
private:
  TString fType; ///< type of this group (e.g. PADS, MANU, PCB, etc...)
  Bool_t fIsResponder; ///< whether we are responding to mouse events
  Bool_t fIsVisible; ///< whether we are visible
  AliMUONVTrackerData* fData; ///< the data we plot (can be 0x0)
  Int_t fDataIndex; ///< the index of the data to plot (can be -1 if data=0x0)
  Double_t fDataMin; ///< min data
  Double_t fDataMax; ///< max data
  TObjArray* fPainters; ///< painters of this group
  Int_t fDepth; ///< depth in the hierarchy of painters
  Bool_t fIsOutlined; ///< whether we should be outlined
  
  ClassDef(AliMUONPainterGroup,1) // Group of AliMUONVPainter
};

#endif
