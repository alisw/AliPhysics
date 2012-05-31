#ifndef ALIMUONSPARSEHISTO_H
#define ALIMUONSPARSEHISTO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONSparseHisto
/// \brief A very memory compact histogram to hold some tracker distributions
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONSparseHisto : public TObject
{
public:
  
  enum 
  {
    kUnderflow = BIT(20),
    kOverflow = BIT(21)
  };
  
  AliMUONSparseHisto(Double_t xmin=0.0, Double_t xmax=4096.0); 
  AliMUONSparseHisto(const AliMUONSparseHisto& rhs);
  AliMUONSparseHisto& operator=(const AliMUONSparseHisto& rhs);
  
  virtual ~AliMUONSparseHisto();
  
  Bool_t Add(const AliMUONSparseHisto& h);
  
  /// Whether this histogram has underflow values 
  /// (no way to know the number of underflow, though)
  Bool_t HasUnderflow() const { return TestBit(kUnderflow); }
  
  /// Whether this histogram has overflow values  
  /// (no way to know the number of underflow, though)
  Bool_t HasOverflow() const { return TestBit(kOverflow); }
  
  Int_t Fill(Double_t value);
  
  /// Return number of bins we hold
  Int_t GetNbins() const { return fNbins; }

  Double_t GetBinCenter(Int_t bin) const;
  
  Int_t GetBinContent(Int_t bin) const;
  
  virtual void Print(Option_t* opt="") const;
  
  virtual void Clear(Option_t* opt="");
    
  Int_t Find(Int_t binCenter) const;
  
  virtual void Copy(TObject& object) const;

  /// Return max value of bincenter
  Double_t Xmax() const { return fXmax; }
  
  /// Return min value of bincenter
  Double_t Xmin() const { return fXmin; }
  
  /// Number of bits used to code the x-value of the histogram
  Int_t Nbits() const { return 12; }
  
private:
  
  UInt_t Encode(Int_t binCenter, Int_t binContent) const;
  
  Double_t DecodeValue(Int_t value) const;
  
  Int_t EncodeValue(Double_t value) const;
  
  UInt_t GetBin(Int_t i) const;
  
  Int_t BinCenter(UInt_t x) const;
  
  Int_t BinContent(UInt_t x) const;
  
  void Expand();
  
  /// Conversion factor to go from float to int value (for bin content)
  Double_t Factor() const { return fFactor; } 
  
private:

  Int_t fNbins;  ///< number of bins we hold

  /// compacted content = (bin,value)
  UInt_t* fArray; //[fNbins] compacted content = (bin,value)

  Double_t fXmin; ///< min value of bincenter
  Double_t fXmax; ///< max value of bincenter
  
  Double_t fFactor; ///< to go from double to int
  
  ClassDef(AliMUONSparseHisto,2) // Sparse histogram-like class for ADC distributions
};

#endif
