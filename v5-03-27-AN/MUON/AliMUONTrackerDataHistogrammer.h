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
  AliMUONTrackerDataHistogrammer(const AliMUONVTrackerData& data, 
                                 Int_t externalDim,
                                 Int_t internalDim=-1);
  virtual ~AliMUONTrackerDataHistogrammer();
  
  static TH1* CreateHisto(const AliMUONVPainter& painter, 
                          Int_t externalDim,
                          Int_t internalDim);

  TH1* CreateChannelHisto(Int_t detElemId, Int_t manuId, 
                          Int_t manuChannel) const;
    
  /// Whether we are working with internal dimensions or external ones.
  Bool_t IsInternalMode() const { return fInternalDim >=0; }
  
private:

  TH1* CreateManuHisto(Int_t detElemId, Int_t manuId, Int_t nbins, Double_t xmin, Double_t xmax) const;
  
  TH1* CreateHisto(const char* basename, Int_t nbins, Double_t xmin, Double_t xmax) const;
  
  void GetDataRange(const TObjArray& manuList, Double_t& xmin, Double_t& xmax) const;
  
  void Add(TH1& h, const AliMUONSparseHisto& sh) const;
  
  void AddBusPatchHisto(TH1& h, Int_t busPatchId) const;
  
  void AddDEHisto(TH1& h, Int_t detElemId) const;
  
  void AddManuHisto(TH1& h, Int_t detElemId, Int_t manuId) const;
  
private:
  const AliMUONVTrackerData& fkData; ///< data we'll histogram 
  Int_t fExternalDim; ///< (external) dimension we'll histogram
  Int_t fInternalDim; ///< (internal) dimension we'll make histogram for
  
  ClassDef(AliMUONTrackerDataHistogrammer,2) // Make histograms from VTrackerData
};

#endif
