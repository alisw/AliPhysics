#ifndef ALIMUONSPARSEHISTO_H
#define ALIMUONSPARSEHISTO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONSparseHisto
/// \brief A very memory compact histogram to hold tracker ADC distributions
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONSparseHisto : public TObject
{
public:
  AliMUONSparseHisto();
  AliMUONSparseHisto(const AliMUONSparseHisto& rhs);
  AliMUONSparseHisto& operator=(const AliMUONSparseHisto& rhs);
  
  virtual ~AliMUONSparseHisto();
  
  Int_t Fill(Int_t adc);
  
  /// Return number of bins we hold
  Int_t GetNbins() const { return fNbins; }
  
  Int_t GetBinContent(Int_t bin) const;
  
  virtual void Print(Option_t* opt="") const;
  
  virtual void Clear(Option_t* opt="");
  
  void Decode(Int_t value, Int_t& adc, Int_t& count) const;
  
  Int_t Encode(Int_t adc, Int_t count) const;
  
  Int_t Find(Int_t adc) const;
  
  virtual void Copy(TObject& object) const;

private:
      
    void Expand();
  
private:

  Int_t fNbins;  ///< number of bins we hold

  /// compacted content = (bin,value)
  Int_t* fArray; //[fNbins] compacted content = (bin,value)
  
  ClassDef(AliMUONSparseHisto,1) // Sparse histogram-like class for ADC distributions
};

#endif
