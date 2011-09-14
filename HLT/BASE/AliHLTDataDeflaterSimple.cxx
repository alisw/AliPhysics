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

/// @file   AliHLTDataDeflaterSimple.cxx
/// @author Matthias Richter
/// @date   2011-08-10
/// @brief  Simple deflater implementation storing frequent values below a
///         maximum value with a reduced bit number and others with the full
///         number of bits.

#include "AliHLTDataDeflaterSimple.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TMath.h"
#include <memory>
#include <algorithm>
#include <iostream>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataDeflaterSimple)

AliHLTDataDeflaterSimple::AliHLTDataDeflaterSimple()
  : AliHLTDataDeflater()
  , fParameterDefinitions()
  , fHistograms(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  if (fHistograms) fHistograms->SetOwner(kTRUE);
}

AliHLTDataDeflaterSimple::~AliHLTDataDeflaterSimple()
{
  // destructor
  Clear();

  if (fHistograms) {
    delete fHistograms;
  }
  fHistograms=NULL;

}

int AliHLTDataDeflaterSimple::AddParameterDefinition(const char* name, int bitLength, int reducedBitLength)
{
  /// add a parameter definition to the configuration, return reference id
  fParameterDefinitions.push_back(AliHLTDataDeflaterParameter(name, bitLength, reducedBitLength));
  if (fHistograms) {
    if (!fHistograms->FindObject(name)) {
      fHistograms->Add(new TH1I(name, name, 100, 0, 0x1<<bitLength));
    }
  }
  return fParameterDefinitions.size()-1;
}

int AliHLTDataDeflaterSimple::AddHistogram(TH1* h)
{
  /// add a histogram for deflater statistic of the corresponding parameter
  if (!fHistograms) fHistograms=new TObjArray;
  if (!fHistograms) return -ENOMEM;
  if (h!=NULL && fHistograms->FindObject(h->GetName())) {
    HLTWarning("parameter with name %s already existing, skipping histogram", h->GetName());
    return -EEXIST;
  }
  if (h) fHistograms->Add(h);
  return 0;
}

bool AliHLTDataDeflaterSimple::OutputParameterBits( int memberId, AliHLTUInt64_t const & value )
{
  // write bit pattern of a member to the current byte and position
  if (memberId>=(int)fParameterDefinitions.size()) return false;

  AliHLTUInt32_t switchBit=fParameterDefinitions[memberId].SwitchBit(value); // 0 -> reduced, 1 -> full
  AliHLTUInt64_t v=fParameterDefinitions[memberId].Value(value);
  AliHLTUInt32_t length=fParameterDefinitions[memberId].ValueLength(value);
  fParameterDefinitions[memberId].IncrementBitCount(value);

  if (fHistograms) {
    TObject* o=fHistograms->FindObject(fParameterDefinitions[memberId].GetName());
    if (o) {
      TH1* h=dynamic_cast<TH1*>(o);
      if (h) {
	h->Fill(v);
      }
    }
  }

  if (!OutputBit(switchBit)) return false;
  return OutputBits(v, length);
}

void AliHLTDataDeflaterSimple::Clear(Option_t * option)
{
  // internal cleanup
  TH2F* hParameterCompression=NULL;
  TH2F* hParameterByteSaving=NULL;
  if (fHistograms) {
    int bins=fParameterDefinitions.size();
    TObject* o=NULL;
    o=fHistograms->FindObject("ParameterCompression");
    if (o) {
      hParameterCompression=dynamic_cast<TH2F*>(o);
    } else {
      hParameterCompression=new TH2F("ParameterCompression", "ParameterCompression", bins, 0, bins, 100, 0., 1.1);
      if (hParameterCompression) fHistograms->Add(hParameterCompression);
    }
    /*
    o=fHistograms->FindObject("ParameterByteSaving");
    if (o) {
      hParameterByteSaving=dynamic_cast<TH2F*>(o);
    } else {
      hParameterByteSaving=new TH2F("ParameterByteSaving", "ParameterByteSaving", bins, 0, bins, 10, 0., 1.1);
      if (hParameterByteSaving) fHistograms->Add(hParameterByteSaving);
    }
    */
  }
  unsigned i=0;
  for (vector<AliHLTDataDeflaterParameter>::iterator m=fParameterDefinitions.begin();
       m!=fParameterDefinitions.end(); m++, i++) {
    int bitLength=m->GetBitLength();
    int valueCount=m->GetValueCount();
    if (bitLength==0 || valueCount==0) continue;
    float ratio=(float)m->GetBitCount();
    ratio/=bitLength*valueCount;
    if (hParameterCompression)
      hParameterCompression->Fill(i, ratio);
    ratio=(1-ratio)*valueCount*bitLength/8;
    if (hParameterByteSaving)
      hParameterByteSaving->Fill(i, ratio);

    m->ResetBitCount();
  }
  AliHLTDataDeflater::Clear(option);
}

void AliHLTDataDeflaterSimple::SaveAs(const char *filename,Option_t */*option*/) const
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
	  ) continue;
      TH1* h=reinterpret_cast<TH1*>(fHistograms->At(i));
      if (!h) continue;
      float entropy=CalcEntropy(h);
      if (entropy<0) continue;
      TString title=h->GetTitle();
      title+=Form(" entropy %.2f", entropy);
      h->SetTitle(title);
    }
    fHistograms->Write();
  }

  file->Close();
}

float AliHLTDataDeflaterSimple::CalcEntropy(TH1* histo, const char* /*option*/, int mode)
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

void AliHLTDataDeflaterSimple::Print(Option_t *option) const
{
  // print info
  Print(cout, option);
}

void AliHLTDataDeflaterSimple::Print(ostream& out, Option_t *option) const
{
  // print to stream
  out << "AliHLTDataDeflaterSimple:" << endl;
  AliHLTUInt64_t bitCount=0;
  AliHLTUInt64_t fullSize=0;
  for (vector<AliHLTDataDeflaterParameter>::const_iterator m=fParameterDefinitions.begin();
       m!=fParameterDefinitions.end(); m++) {
    cout << "   "; m->Print(option);
    bitCount+=m->GetBitCount();
    fullSize+=m->GetValueCount()*m->GetBitLength();
  }
  out << " total: " << bitCount << "/" << fullSize << " " << (fullSize>0?float(bitCount)/fullSize:0.0) << endl;
}

void AliHLTDataDeflaterSimple::AliHLTDataDeflaterParameter::Print(const char* /*option*/) const
{
  // print info
  cout << fName << " (" << fFullBitLength << "," << fReducedBitLength << "): "
       << fValueCount << " entries  "
       << fBitCount << "/" << fFullBitLength*fValueCount;
  if (fFullBitLength && fValueCount) {
    cout << " " << float(fBitCount)/(fValueCount*fFullBitLength);
  }
  cout << endl;
}

ostream& operator<<(ostream &out, const AliHLTDataDeflaterSimple& me)
{
  me.Print(out);
  return out;
}
