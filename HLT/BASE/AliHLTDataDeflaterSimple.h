//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTDATADEFLATERSIMPLE_H
#define ALIHLTDATADEFLATERSIMPLE_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTDataDeflaterSimple.h
/// @author Matthias Richter
/// @date   2011-08-10
/// @brief  Data deflater class storing only necessary bits
/// @note   Code original from AliHLTTPCCompModelDeflater

#include "AliHLTDataDeflater.h"
#include <vector>
#include <string>

class TObjArray;
class TH1;

/**
 * @class AliHLTDataDeflaterSimple
 * Simple deflater implementation storing frequent values below a
 * maximum value with a reduced bit number and others with the full
 * number of bits. The reduced value is indicated by a preceeding '0'
 * and the full bit by '1'. The algorithm can be applied to data with an
 * occurrence distribution peaking close to zero and having less frequent
 * occurrence at higher values.
 *
 * @ingroup alihlt_base
 */
class AliHLTDataDeflaterSimple : public AliHLTDataDeflater
{
public:
  /// standard constructor
  AliHLTDataDeflaterSimple();
  /// destructor
  ~AliHLTDataDeflaterSimple();

  /// @class AliHLTDataDeflaterParameter definition of parameters
  class AliHLTDataDeflaterParameter {
  public:
    AliHLTDataDeflaterParameter()
      : fName(), fFullBitLength(0), fReducedBitLength(0)
      , fMax(0), fMaxReduced(0)
      , fMask(0), fMaskReduced(0)
      , fValueCount(0), fBitCount(0) {}

    AliHLTDataDeflaterParameter(const char* name, int length, int reduced)
      : fName(name), fFullBitLength(length), fReducedBitLength(reduced)
      , fMax((((AliHLTUInt64_t)0x1)<<length)-1), fMaxReduced((((AliHLTUInt64_t)0x1)<<reduced)-1)
      , fMask(fMax), fMaskReduced(fMaxReduced) 
      , fValueCount(0), fBitCount(0) {}

    AliHLTDataDeflaterParameter(const AliHLTDataDeflaterParameter& src)
      : fName(src.fName), fFullBitLength(src.fFullBitLength), fReducedBitLength(src.fReducedBitLength)
      , fMax(src.fMax), fMaxReduced(src.fMaxReduced)
      , fMask(src.fMask), fMaskReduced(src.fMaskReduced) 
      , fValueCount(0), fBitCount(0) {}

    AliHLTDataDeflaterParameter& operator=(const AliHLTDataDeflaterParameter& src) {
      fName=src.fName; fFullBitLength=src.fFullBitLength; fReducedBitLength=src.fReducedBitLength;
      fMax=src.fMax; fMaxReduced=src.fMaxReduced;
      fMask=src.fMask; fMaskReduced=src.fMaskReduced;
      fValueCount=src.fValueCount; fBitCount=src.fBitCount;
      return *this;
    }

    ~AliHLTDataDeflaterParameter() {}

    const char* GetName() const {return fName.c_str();}
    AliHLTUInt64_t Value(const AliHLTUInt64_t& value) const {
      return value>fMax?fMax:value;
    }
    AliHLTUInt32_t ValueLength(const AliHLTUInt64_t& value) const{
      return value>fMaxReduced?fFullBitLength:fReducedBitLength;
    }
    AliHLTUInt32_t SwitchBit(const AliHLTUInt64_t& value) const {
      return value>fMaxReduced;
    }
    const int& GetBitLength() const {return fFullBitLength;}
    const int& GetReducedBitLength() const {return fReducedBitLength;}
    const AliHLTUInt64_t& GetMax() const {return fMax;}
    const AliHLTUInt64_t& GetMaxReduced() const {return fMaxReduced;}
    const AliHLTUInt64_t& GetMask() const {return fMask;}
    const AliHLTUInt64_t& GetReducedMask() const {return fMaskReduced;}

    const AliHLTUInt32_t& GetBitCount() const {return fBitCount;}
    const AliHLTUInt32_t& GetValueCount() const {return fValueCount;}
    void IncrementBitCount(const AliHLTUInt64_t& value) {
      fBitCount+=(value>fMaxReduced?fFullBitLength:fReducedBitLength)+1;
      fValueCount++;
    }
    void ResetBitCount() {fValueCount=0; fBitCount=0;}
    void Print(const char* option="") const;

  private:
    std::string fName; //!
    int fFullBitLength; //!
    int fReducedBitLength; //!
    AliHLTUInt64_t fMax; //!
    AliHLTUInt64_t fMaxReduced; //!
    AliHLTUInt64_t fMask; //!
    AliHLTUInt64_t fMaskReduced; //!
    AliHLTUInt32_t fValueCount; //!
    AliHLTUInt32_t fBitCount; //!
  };

  /// add a parameter definition to the configuration, return reference id
  int AddParameterDefinition(const char* name, int bitLength, int reducedBitLength);

  /// add a histogram for deflater statistic of the corresponding parameter
  int AddHistogram(TH1* h);

  /// inherited from AliHLTDataDeflater: write bit pattern according to configuration
  virtual bool OutputParameterBits( int parameterId, AliHLTUInt64_t const & value );

  /// clear the object and reset pointer references
  virtual void Clear(Option_t * /*option*/ ="");

  /// print info
  virtual void Print(Option_t *option="") const;

  /// print info
  virtual void Print(ostream& out, Option_t *option="") const;

  /// safe statistics histograms to file
  virtual void SaveAs(const char *filename="",Option_t *option="") const;

  /// DataDeflaterSimple has deflater version 1
  virtual int GetDeflaterVersion() const {return 1;}

 protected:
 private:
  /// copy constructor prohibited
  AliHLTDataDeflaterSimple(const AliHLTDataDeflaterSimple&);
  /// assignment operator prohibited
  AliHLTDataDeflaterSimple& operator=(const AliHLTDataDeflaterSimple&);

  vector<AliHLTDataDeflaterParameter> fParameterDefinitions; //!

  TObjArray* fHistograms; //! list of histograms for parameters

  ClassDef(AliHLTDataDeflaterSimple, 0)
};

ostream& operator<<(ostream &out, const AliHLTDataDeflaterSimple& me);

#endif
