// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTDataDeflater.cxx
/// @author Matthias Richter, Timm Steinbeck
/// @date   2011-08-10
/// @brief  Data deflater class storing only necessary bits
/// @note   Code original from AliHLTTPCCompModelDeflater

#include "AliHLTDataDeflater.h"
#include "AliHLTErrorGuard.h"
#include "TObjArray.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include <memory>
#include <algorithm>
#include <iostream>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataDeflater)

AliHLTDataDeflater::AliHLTDataDeflater()
  : AliHLTLogging()
  , fBitDataCurrentWord(0)
  , fBitDataCurrentPosInWord(0)
  , fBitDataCurrentOutput(NULL)
  , fBitDataCurrentOutputStart(NULL)
  , fBitDataCurrentOutputEnd(NULL)
  , fHistograms(NULL)
  , fParameterCompression(NULL)
  , fParameterSize(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  if (fHistograms) fHistograms->SetOwner(kTRUE);
}

AliHLTDataDeflater::~AliHLTDataDeflater()
{
  // destructor
  Clear();

  if (fHistograms)
    delete fHistograms;
  fHistograms=NULL;
  if (fParameterCompression)
    delete fParameterCompression;
  fParameterCompression=NULL;
  if (fParameterSize)
    delete fParameterSize;
  fParameterSize=NULL;
}

int AliHLTDataDeflater::InitBitDataOutput( AliHLTUInt8_t* output, UInt_t outputSize)
{
  // init the target buffer
  fBitDataCurrentWord = 0;
  fBitDataCurrentPosInWord = 7;
  fBitDataCurrentOutput = fBitDataCurrentOutputStart = output;
  fBitDataCurrentOutputEnd = output+outputSize;

  return 0;
}

void AliHLTDataDeflater::CloseBitDataOutput()
{
  // pad to full byte and clear internal pointer references
  Pad8Bits();
  fBitDataCurrentWord=0;
  fBitDataCurrentPosInWord=0;
  fBitDataCurrentOutput=NULL;
  fBitDataCurrentOutputStart=NULL;
  fBitDataCurrentOutputEnd=NULL;
}

AliHLTUInt8_t AliHLTDataDeflater::GetCurrentOutputByte( Int_t offset ) const
{
  // get the current byte
  if ( !offset )
    return fBitDataCurrentWord;
  else
    return *(fBitDataCurrentOutput+offset);
}

bool AliHLTDataDeflater::OutputBit( AliHLTUInt32_t const & value )
{
  // write one bit to the current byte and position
  if ( fBitDataCurrentOutput>=fBitDataCurrentOutputEnd )
    return false;
  fBitDataCurrentWord |= (value & 1) << fBitDataCurrentPosInWord;
  if ( fBitDataCurrentPosInWord )
    fBitDataCurrentPosInWord--;
  else {
    *fBitDataCurrentOutput = fBitDataCurrentWord;
    fBitDataCurrentPosInWord = 7;
    fBitDataCurrentOutput++;
    fBitDataCurrentWord = 0;
  }
  return true;
}

bool AliHLTDataDeflater::OutputBits( AliHLTUInt64_t const & value, UInt_t const & bitCount )
{
  // write bit pattern to the current byte and position
  if ( bitCount>64 ) {
    HLTFatal( "Internal error: Attempt to write more than 64 bits (%u)", (unsigned)bitCount );
    return false;
  }
  UInt_t bitsToWrite=bitCount;
  UInt_t curBitCount;
  while ( bitsToWrite>0 ) {
    if ( fBitDataCurrentOutput>=fBitDataCurrentOutputEnd )
      return false;
#if 1
    if ( bitsToWrite >= fBitDataCurrentPosInWord+1 )
      curBitCount = fBitDataCurrentPosInWord+1;
    else
      curBitCount = bitsToWrite;
    fBitDataCurrentWord |= ( (value >> (bitsToWrite-curBitCount)) & ((1<<curBitCount)-1) ) << (fBitDataCurrentPosInWord+1-curBitCount);
    if ( fBitDataCurrentPosInWord < curBitCount )
      {
	*fBitDataCurrentOutput = fBitDataCurrentWord;
	fBitDataCurrentPosInWord = 7;
	fBitDataCurrentOutput++;
	fBitDataCurrentWord = 0;
      }
    else
      fBitDataCurrentPosInWord -= curBitCount;
    bitsToWrite -= curBitCount;

#else
    AliHLTUInt8_t curValue;
    if ( bitsToWrite>=8 )
      {
	curBitCount=8;
	curValue = (value >> bitsToWrite-8) & 0xFF;
	bitsToWrite -= 8;
      }
    else
      {
	curBitCount=bitsToWrite;
	curValue = value & ( (1<<bitsToWrite)-1 );
	bitsToWrite = 0;
      }
    if ( fBitDataCurrentPosInWord+1>curBitCount )
      {
	fBitDataCurrentWord |= curValue << (fBitDataCurrentPosInWord-curBitCount+1);
	fBitDataCurrentPosInWord -= curBitCount;
      }
    else if ( fBitDataCurrentPosInWord+1==curBitCount )
      {
	fBitDataCurrentWord |= curValue;
	*fBitDataCurrentOutput = fBitDataCurrentWord;
	fBitDataCurrentPosInWord = 7;
	fBitDataCurrentOutput++;
	fBitDataCurrentWord = 0;
      }
    else
      {
	const UInt_t first = fBitDataCurrentPosInWord+1; // Number of bits for first block
	const UInt_t second = curBitCount-first; // Number of bits for second block
	fBitDataCurrentWord |= ( curValue >> second ) & ((1<<first)-1);
	*fBitDataCurrentOutput = fBitDataCurrentWord;
	fBitDataCurrentOutput++;
	if ( fBitDataCurrentOutput>=fBitDataCurrentOutputEnd )
	  return false;
	fBitDataCurrentWord = curValue & ((1<<second)-1) << (8-second);
	fBitDataCurrentPosInWord = 7-second;
      }
#endif
  }
  return true;
}

bool AliHLTDataDeflater::OutputBits( std::bitset<64> const & value, UInt_t const & bitCount )
{
  // write bit pattern to the current byte and position
  if ( bitCount>64 ) {
    HLTFatal( "Internal error: Attempt to write more than 64 bits (%u)", (unsigned)bitCount );
    return false;
  }
  static const std::bitset<64> mask8bit(255ul);
  UInt_t bitsToWrite=bitCount;
  UInt_t curBitCount;
  while ( bitsToWrite>0 ) {
    if ( fBitDataCurrentOutput>=fBitDataCurrentOutputEnd )
      return false;
    if ( bitsToWrite >= fBitDataCurrentPosInWord+1 )
      curBitCount = fBitDataCurrentPosInWord+1;
    else
      curBitCount = bitsToWrite;
    std::bitset<64> valwrite=(value >> (bitsToWrite-curBitCount)) & mask8bit;
    fBitDataCurrentWord |= ( valwrite.to_ulong() & ((1<<curBitCount)-1) ) << (fBitDataCurrentPosInWord+1-curBitCount);
    if ( fBitDataCurrentPosInWord < curBitCount )
      {
	*fBitDataCurrentOutput = fBitDataCurrentWord;
	fBitDataCurrentPosInWord = 7;
	fBitDataCurrentOutput++;
	fBitDataCurrentWord = 0;
      }
    else
      fBitDataCurrentPosInWord -= curBitCount;
    bitsToWrite -= curBitCount;
  }
  return true;
}

void AliHLTDataDeflater::Pad8Bits()
{
  // finish the current word
  if ( fBitDataCurrentPosInWord==7 )
    return;
  *fBitDataCurrentOutput = fBitDataCurrentWord;
  fBitDataCurrentPosInWord = 7;
  fBitDataCurrentOutput++;
  fBitDataCurrentWord = 0;
}

bool AliHLTDataDeflater::OutputBytes( AliHLTUInt8_t const * data, UInt_t const & byteCount )
{
  // write sequence of bytes
  Pad8Bits();
  if ( fBitDataCurrentOutput+byteCount>fBitDataCurrentOutputEnd )
    return false;
  memcpy( fBitDataCurrentOutput, data, byteCount );
  fBitDataCurrentOutput += byteCount;
  return true;
}

bool AliHLTDataDeflater::OutputParameterBits( int parameterId, AliHLTUInt64_t const & value )
{
  // write bit pattern of a member to the current byte and position
  return OutputParameterBits(parameterId, value, 0);
}

bool AliHLTDataDeflater::OutputParameterBits( int /*(parameterId*/, AliHLTUInt64_t const & /*value*/ , int /*lengthOffset*/)
{
  // write bit pattern of a member to the current byte and position
  ALIHLTERRORGUARD(1,"method needs to be implemented in child class");
  return false;
}

void AliHLTDataDeflater::Clear(Option_t * /*option*/)
{
  // internal cleanup
}

void AliHLTDataDeflater::Print(Option_t *option) const
{
  // print info
  Print(cout, option);
}

void AliHLTDataDeflater::Print(ostream& out, Option_t */*option*/) const
{
  // print to stream
  out << "AliHLTDataDeflater: " << endl;
}

ostream& operator<<(ostream &out, const AliHLTDataDeflater& me)
{
  me.Print(out);
  return out;
}

int AliHLTDataDeflater::EnableStatistics()
{
  /// enable statistics accounting
  if (!fHistograms) {
    fHistograms=new TObjArray;
    if (!fHistograms) return -ENOMEM;
    fHistograms->SetOwner(kTRUE);
  }
  return 0;
}

int AliHLTDataDeflater::AddHistogram(int id, const char* name, int bitWidth, TH1* h)
{
  /// add a histogram for deflater statistic of the corresponding parameter
  /// a histogram is created with default settings if h is not specified; if
  /// provided, the ownership goes over to the base class
  if (!fHistograms) {
    fHistograms=new TObjArray;
    if (!fHistograms) return -ENOMEM;
    fHistograms->SetOwner(kTRUE);
  }
  if (id>=0 && fHistograms->GetEntriesFast()>id && fHistograms->At(id)!=NULL) {
    HLTWarning("parameter with id %d has existing object (%s), skipping histogram %s",
	       id, fHistograms->At(id)->GetName(), h?h->GetName():name);
    return -EEXIST;
  }
  if (id<0 && h!=NULL && fHistograms->FindObject(h->GetName())) {
    HLTWarning("parameter with name %s already existing, skipping histogram", h->GetName());
    return -EEXIST;
  }
  if (!h)
    h=new TH1I(name, name, 100, 0, 1<<bitWidth);
  if (!h) return -ENOMEM;
  fHistograms->Add(h);
  return 0;
}

int AliHLTDataDeflater::FillStatistics(int id, AliHLTUInt64_t value, unsigned length, float codingWeight)
{
  /// fill statistics for a parameter
  if (!fHistograms ||
      fHistograms->GetEntriesFast()<=id ||
      id<0) return 0;

  if (value<(~(AliHLTUInt64_t)0)) {
  TObject* o=fHistograms->At(id);
  if (o) {
    TH1* h=dynamic_cast<TH1*>(o);
    if (h) {
      h->Fill(value);
    }
  }
  }

  if (!fParameterCompression) {
    int bins=fHistograms->GetEntriesFast();
    fParameterCompression=new TH2D("ParameterCompression", "ParameterCompression", bins, 0, bins, 1000, 0., 5.0);
  }
  if (fParameterCompression && codingWeight>=.0) {
    fParameterCompression->Fill(id, codingWeight);
  }
  if (!fParameterSize) {
    int bins=fHistograms->GetEntriesFast();
    fParameterSize=new TH2D("ParameterSize", "ParameterSize", bins, 0, bins, 1000, 0., 64.0);
  }
  if (fParameterSize && length>0) {
    fParameterSize->Fill(id, length);
  }

  return 0;
}

void AliHLTDataDeflater::SaveAs(const char *filename,Option_t */*option*/) const
{
  // safe histograms to file
  std::auto_ptr<TFile> file(TFile::Open(filename, "RECREATE"));
  if (!file.get() || file->IsZombie()) {
    HLTError("can not open file %s", filename);;
    return;
  }
  file->cd();
  if (fHistograms) {
    for (int i=0; i<fHistograms->GetEntries(); i++) {
      if (fHistograms->At(i)==NULL || 
	  !fHistograms->At(i)->InheritsFrom("TH1") ||
	  fHistograms->At(i)->InheritsFrom("TH2") ||
	  fHistograms->At(i)->InheritsFrom("TH3")
	  ) continue; // only TH1 objects in the list
      TH1* h=reinterpret_cast<TH1*>(fHistograms->At(i));
      if (!h) continue;
      float entropy=CalcEntropy(h);
      if (entropy<0) continue;
      TString title=h->GetTitle();
      title+=Form(" entropy %.2f", entropy);
      h->SetTitle(title);
    }
    fHistograms->Write();
    if (fParameterCompression)
      fParameterCompression->Write();
    if (fParameterSize)
      fParameterSize->Write();
  }

  file->Close();
}

float AliHLTDataDeflater::CalcEntropy(TH1* histo, const char* /*option*/, int mode)
{

  if (!histo) return -1000.;

  float l2=TMath::Log(2.0);
  float integral=histo->Integral(0,histo->GetNbinsX());
  int centerbin=mode*histo->GetNbinsX()/2;
  int nofBins=histo->GetNbinsX()-centerbin;
  float entropy=0.0;
  for (int offset=0; offset<nofBins; offset++) {
    float abundance=histo->GetBinContent(offset);
    if (abundance<1.0) continue;
    entropy += (- (Double_t) abundance / (Double_t) integral ) * log( ( (Double_t) abundance / (Double_t) integral )) / (l2);
  }

  return entropy;
}
