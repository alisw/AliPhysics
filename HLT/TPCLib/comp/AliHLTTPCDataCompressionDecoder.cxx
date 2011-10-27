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

/// @file   AliHLTTPCDataCompressionDecoder.cxx
/// @author Matthias Richter
/// @date   2011-10-04
/// @brief  Generic decoder class for compressed TPC data, works on a container
///         class implementation which fills the actual target data struct

#include "AliHLTTPCDataCompressionDecoder.h"
#include "AliHLTDataInflaterSimple.h"
#include "AliHLTDataInflaterHuffman.h"
#include "TList.h"
#include <memory>

ClassImp(AliHLTTPCDataCompressionDecoder)

AliHLTTPCDataCompressionDecoder::AliHLTTPCDataCompressionDecoder()
  : fVerbosity(0)
{
  /// constructor
}

AliHLTTPCDataCompressionDecoder::~AliHLTTPCDataCompressionDecoder()
{
  ///destructor
}

AliHLTDataInflater* AliHLTTPCDataCompressionDecoder::CreateInflater(int deflater, int mode) const
{
  // create the inflater for the specified mode
  vector<AliHLTTPCDefinitions::AliClusterParameterId_t> parameterids;
  switch (mode) {
  case 1:
    parameterids.push_back(AliHLTTPCDefinitions::kPadRow );
    parameterids.push_back(AliHLTTPCDefinitions::kPad    );
    parameterids.push_back(AliHLTTPCDefinitions::kTime   );
    parameterids.push_back(AliHLTTPCDefinitions::kSigmaY2);
    parameterids.push_back(AliHLTTPCDefinitions::kSigmaZ2);
    parameterids.push_back(AliHLTTPCDefinitions::kCharge );
    parameterids.push_back(AliHLTTPCDefinitions::kQMax   );
    break;
  case 2:
    parameterids.push_back(AliHLTTPCDefinitions::kResidualPad );
    parameterids.push_back(AliHLTTPCDefinitions::kResidualTime);
    parameterids.push_back(AliHLTTPCDefinitions::kSigmaY2);
    parameterids.push_back(AliHLTTPCDefinitions::kSigmaZ2);
    parameterids.push_back(AliHLTTPCDefinitions::kCharge );
    parameterids.push_back(AliHLTTPCDefinitions::kQMax   );
    break;
  default:
    HLTError("invalid mode %d for inflater initialization", mode);
  }

  switch (deflater) {
  case 1:
    {
      std::auto_ptr<AliHLTDataInflaterSimple> inflatersimple(new AliHLTDataInflaterSimple);
      if (!inflatersimple.get()) return NULL;
      for (vector<AliHLTTPCDefinitions::AliClusterParameterId_t>::const_iterator id=parameterids.begin();
	   id!=parameterids.end(); id++) {
	const AliHLTTPCDefinitions::AliClusterParameter& parameter=AliHLTTPCDefinitions::fgkClusterParameterDefinitions[*id];
	if (inflatersimple->AddParameterDefinition(parameter.fName,
						   parameter.fBitLength,
						   parameter.fOptional)<0) {
	  HLTError("error adding parameter definition %s to inflater", parameter.fName);
	  return NULL;
	}
      }
      return inflatersimple.release();
    }
    break;
  case 2:
    {
      std::auto_ptr<AliHLTDataInflaterHuffman> inflaterhuffman(new AliHLTDataInflaterHuffman);
      if (!inflaterhuffman.get()) return NULL;
      TString cdbPath("HLT/ConfigTPC/TPCDataCompressorHuffmanTables");
      TObject* pConf=AliHLTMisc::Instance().ExtractObject(AliHLTMisc::Instance().LoadOCDBEntry(cdbPath));
      if (!pConf) {
	HLTError("can not load configuration object %s", cdbPath.Data());
	return NULL;
      }
      if (dynamic_cast<TList*>(pConf)==NULL) {
	HLTError("huffman table configuration object of inconsistent type");
	return NULL;
      }
      inflaterhuffman->InitDecoders(dynamic_cast<TList*>(pConf));
      for (vector<AliHLTTPCDefinitions::AliClusterParameterId_t>::const_iterator id=parameterids.begin();
	   id!=parameterids.end(); id++) {
	const AliHLTTPCDefinitions::AliClusterParameter& parameter=AliHLTTPCDefinitions::fgkClusterParameterDefinitions[*id];
	if (inflaterhuffman->AddParameterDefinition(parameter.fName,
						    parameter.fBitLength)<0) {
	  HLTError("error adding parameter definition %s to inflater", parameter.fName);
	  return NULL;
	}
      }
      return inflaterhuffman.release();
    }
    break;
  default:
    HLTError("unknown inflater requested %d", deflater);
  }
  return NULL;
}
