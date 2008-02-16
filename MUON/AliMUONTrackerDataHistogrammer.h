#ifndef ALIMUONTRACKERDATAHISTOGRAMMER_H
#define ALIMUONTRACKERDATAHISTOGRAMMER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONTrackerDataHistogrammer
/// \brief Make histograms from VTrackerData and VPainter objects.
/// 
// author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVPainter;
class AliMUONVTrackerData;
class AliMUONSparseHisto;
class TH1;

class AliMUONTrackerDataHistogrammer : public TObject
{
public:
  AliMUONTrackerDataHistogrammer(const AliMUONVTrackerData& data, Int_t dim);
  virtual ~AliMUONTrackerDataHistogrammer();
  
  static TH1* CreateHisto(const AliMUONVPainter& painter);
  
  TH1* CreateBusPatchHisto(Int_t busPatchId) const;
  
  TH1* CreateChamberHisto(Int_t chamberId) const;
  
  TH1* CreateChannelHisto(Int_t detElemId, Int_t manuId, 
                          Int_t manuChannel) const;
  
  TH1* CreateDEHisto(Int_t detElemId) const;
  
  TH1* CreateManuHisto(Int_t detElemId, Int_t manuId) const;
  
private:

  TH1* CreateHisto(const char* name) const;
  
  void Add(TH1& h, const AliMUONSparseHisto& sh) const;
  
  void AddBusPatchHisto(TH1& h, Int_t busPatchId) const;
  
  void AddDEHisto(TH1& h, Int_t detElemId) const;
  
  void AddManuHisto(TH1& h, Int_t detElemId, Int_t manuId) const;
  
private:
  const AliMUONVTrackerData& fData; ///< data we'll histogram 
  Int_t fDim; ///< dimension we'll histogram
  
  ClassDef(AliMUONTrackerDataHistogrammer,1) // Make histograms from VTrackerData
};

#endif
